#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script generates the specified grammatical relation tuples
from a plain text corpus using the Stanford parser. The input file is assumed
to be a plain text with one sentence per line. All xml tags will be ingnored.

USAGE: generate_rels.py NUMBER_OF_CORES INPUT_FILE
    NUMBER_OF_CORES       Number of java workers. This number should not be grater
                          than number of CPU cores.
    INPUT_FILE            Path to the input file.
    -h, --help            Print this help.
    -m, --method m        Use method of selecting frame realizations m. The default
                          method is "allrels".
    -s, --starting_line s Start the parsing from line #s. This serves for recovering reasons.
                          When this option is set, the output data is appended to the output 
                          file instead of replacing it.
    -o, --output_file o   Use output file o. Default is the same as the input with ".rel" 
                          suffix.
    -n, --corenlp n       Path to CoreNLP directory. Default is 
                          "../../3rdparty/stanford-corenlp/".
    -x, --max_memory x    Maximum amount of memory for each core in gigabytes. Default is 3. 
"""


import sys
import re
import time
import signal
import getopt
import methods
import threading
import subprocess
import logging
import logger
from Queue import Queue
from logger import DBG
from logger import WARN
from logger import INFO

logger.setLogFile("generate_rels.log")
logger.setLevel(logging.INFO)

class Worker:
    def __init__(self, args):
        self.args = args
        self.server = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            preexec_fn = self.preexec_function)
        self.stdin = self.server.stdin
        self.stdout = self.server.stdout
        DBG("Worker initialized.")

    def preexec_function(self):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    
    def parse(self, inputsentence, recover=False):
        DBG("Parsing sentence: %s" % inputsentence)
        
        try:
            self.stdin.write(inputsentence + "\n")
        except IOError:
            WARN("Error while writing the sentence to CoreNLP.")
            return False
        state = 0
        result = {"deps":[], "lemmas":[]}

        DBG("Worker wrote sentence")

        while(True):
            #DBG("Waiting for line.")
            line = self.stdout.readline()
            #DBG("Read line: '%s'." % line)
            if line == "":
                WARN("Read empty line. CoreNLP server is probably down.")
                return False
            if line.startswith("NLP"):
                state = 1
                continue
            if state == 1:
                result["text"] = line.strip()
                state = 2
                continue
            if line.startswith("[Text="):                
                exp = re.compile('\[([^\]]+)\]')
                matches  = exp.findall(line)
                for s in matches:
                    av = re.split("=| ", s)
                    word = av[1] 
                    lemma = ""
                    attributes = {}
                    for a,v in zip(*[av[2:][x::2] for x in (0, 1)]):
                        if a == "Lemma":
                            lemma = v
                            break   
                    result['lemmas'].append(lemma)
            if line.startswith("(ROOT"):
                state = 3
                continue
            if line.strip() == "" and state == 3:
                state = 4
                continue

            if line.strip() != "" and state == 4:
                line = line.rstrip()
                split_entry = re.split("\(|, ", line[:-1])
                if len(split_entry) == 3:
                    rel, left, right = map(lambda x: x, split_entry)
                    result['deps'].append([rel,left,right])
    
            if line.strip() == "" and state == 4:
                DBG("Worker returned processed sentence.")
                return [result]
        
    def stop(self):
        self.server.kill()
        self.server.wait()

    def __del__(self):
        self.stop()

class Client(threading.Thread):
    def __init__(self, sentence, server):
        self.sentence = sentence
        self.server = server
        self.result = False
        threading.Thread.__init__(self)

    def getSentence(self):
        return self.sentence

    def getResult(self):
        return self.result

    def getServer(self):
        return self.server

    def run(self):
        self.result = self.server.parse(self.sentence)


class Generator:
    def __init__(self):
        self.servers = Queue()
        self.processedLines = 0
        self.processedSentences = 0
        self.threads = Queue()
        self.finished = False
        self.exit = False

    def createWorkers(self, cores, args):
        self.cores = cores
        DBG("Creating %d workers." % cores)
        for s in range(self.cores):
            self.servers.put(Worker(args))

    def producer(self, inputFileName, startingLine):
        INFO("Starting producer.")
        inputFile = open(inputFileName, "r")
        for line in inputFile.xreadlines():
            if self.exit == True:
                break
            self.processedLines += 1
            if self.processedLines < startingLine:
                continue
            line = line.strip()
            if line.startswith("<") or len(line) < 5:
                continue

            server = self.servers.get(True)
            thread = Client(line, server)
            thread.start()
            self.threads.put(thread, True)
            if self.processedLines % (self.cores*10) == 0:
                DBG("Processed lines: %d." % self.processedLines)
                sys.stderr.write("Processed lines: %d.\n" % self.processedLines)
            self.processedSentences += 1

        inputFile.close()
        self.finished = True
        INFO("Exiting producer.")


    def consumer(self, outputFileName, startingLine, method):
        INFO("Starting consumer.")
        outputFile = None
        writtenSentences = 0

        if startingLine == 0:
            outputFile = open(outputFileName, "w")
            outputFile.write(method(None, desc = True) + "\n")
        else:
            outputFile = open(outputFileName, "a")


        while not self.exit and (not self.finished or \
                writtenSentences < self.processedSentences):

            thread = self.threads.get(True)
            DBG("Consumer got worker.")
            thread.join(300) #max. time to process one sentence is five minutes
            server = thread.getServer()
            if (thread.isAlive()):
                WARN("Server timed out. Skipping line and restarting server. " +\
                     "(Sentence: '%s')"  %  thread.getSentence())
                server.server.kill()
                self.servers.put(Worker(server.args), True)
            else:
                DBG("Worker finished.")
                result = thread.getResult()
                tuples = method(result)

                if result == False:
                    WARN("Internal error occured. Restarting server.")
                    server.server.kill()
                    self.servers.put(Worker(server.args), True)
                else:
                    outputFile.write(tuples)
                    self.servers.put(server, True)
                    writtenSentences += 1

        outputFile.close()
        INFO("Exiting consumer.")

    def signalHandler(self, signum, frame):
        INFO("Last processed line: %d." % self.processedLines)
        self.exit = True


    def run(self, inputFileName, outputFileName, method, startingLine = 0):
        
        self.cons_thread = threading.Thread(target=self.consumer,
                args=(outputFileName, startingLine, method))
        self.cons_thread.start()

        self.producer(inputFileName, startingLine)

        self.cons_thread.join()


if __name__ == "__main__":

    cores = 0
    input_file = ""
    starting_line = 0
    output_file = ""
    corenlp = "../../3rdparty/stanford-corenlp/"
    max_memory = 3
    method = methods.parserstruct
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:s:d:c:x:", ["help", "method",
            "starting_line", "output_dir", "corenlp", "max_memory"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        if o in ("-m", "--method"):
            try:
                method = getattr(methods, a)
            except AttributeError, msg:
                print msg
                print "for help use --help"
                sys.exit(2)
        if o in ("-s", "--starting_line"):
            starting_line = int(a)
        if o in ("-o", "--output_file"):
            output_file = a
        if o in ("-c", "--corenlp"):
            corenlp = a
        if o in ("-x", "--max_memory"):
            max_memory = int(a)


    if len(args) != 2:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)
    cores = int(args[0])
    input_file = args[1]

    if output_file == "": output_file = input_file + ".rel"

    if not corenlp.endswith("/"):
        corenlp += "/"    

    jars = ":".join([
        corenlp + "stanford-corenlp-2012-01-08.jar",
        corenlp + "stanford-corenlp-2011-12-27-models.jar",
        corenlp + "xom.jar",
        corenlp + "joda-time.jar"
    ])

    args = ["/usr/bin/java",
            "-cp", 
            jars, 
            "-Xmx%dg" % max_memory, 
            "edu.stanford.nlp.pipeline.StanfordCoreNLP",
            "-annotators",  "tokenize,ssplit,pos,lemma,parse"]

    generator = Generator()
    signal.signal(signal.SIGINT, generator.signalHandler)
    generator.createWorkers(cores, args)
    generator.run(input_file, output_file, method, starting_line)
    sys.exit(0)

