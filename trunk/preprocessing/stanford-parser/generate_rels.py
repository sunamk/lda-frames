#!/usr/bin/python

"""This script generates the specified grammatical relation tuples
from a plain text corpus using the Stanford parser. The input file is assumed
to be a plain text with one sentence per line. All xml tags will be ingnored.

USAGE: generate_rels.py NUMBER_OF_CORES STARTING_PORT_NUMBER INPUT_FILE
    -h, --help            Print this help.
    -m, --method m        Use method of selecting frame realizations m. The default
                          method is "allrels".
    -s, --starting_line s Start the parsing from line #s. This serves for recovering reasons.
                          When this option is set, the output data is appended to the output 
                          file instead of replacing it.
    -o, --output_file o   Use output file o. Default is the same as the input with ".rel" 
                          suffix.
    -n, --host_name n     Connect to host n. Defaut is "127.0.0.1".
"""


import sys
import signal
import getopt
import jsonrpc
import methods
import threading
import time
from Queue import Queue
from simplejson import loads

class Client(threading.Thread):
    def __init__(self, sentence, server):
        self.sentence = unicode(sentence, "utf-8")
        self.server = server
        self.result = False
        threading.Thread.__init__(self)

    def getResult(self):
        return self.result

    def getServer(self):
        return self.server

    def run(self):
        try:
            self.result = loads(self.server.parse(self.sentence, False))
            if self.result == []:
                sys.stderr.write("Empty result. Waiting 10 seconds.\n")
                time.sleep(10)
                self.run()
    
        except jsonrpc.RPCTransportError, msg:
            if str(msg) == "timed out":
                sys.stderr.write(str(msg) + ". Ignoring line.\n")
                self.result = None
            else:
                sys.stderr.write("Transport error. " + str(msg) + "\n")
                
        except jsonrpc.RPCInternalError, msg:
            sys.stderr.write(str(msg))
            self.result = False

        except jsonrpc.RPCParseError, msg:
            sys.stderr.write(str(msg) + ". Ignoring line.\n")
            self.result = None
        


class Generator:
    def __init__(self, cores, starting_port, host_name):
        self.starting_port = starting_port
        self.cores = cores
        self.host_name = host_name
        self.servers = Queue()
        self.processedLines = 0
        self.processedSentences = 0
        self.finished = False
        self.threads = Queue()
        self.exit = False

        for p in range(self.starting_port, self.starting_port + self.cores):
            self.servers.put(jsonrpc.ServerProxy(jsonrpc.JsonRpc20(),
                jsonrpc.TransportTcpIp(addr=(self.host_name, p))))

    def producer(self, inputFileName, startingLine):
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
                sys.stderr.write("Processed lines: %d.\n" % self.processedLines)
            self.processedSentences += 1
        
        inputFile.close()
        self.finished = True

    def consumer(self, outputFileName, startingLine, method):
        outputFile = None
        writtenSentences = 0

        if startingLine == 0:
            outputFile = open(outputFileName, "w")
        else:
            outputFile = open(outputFileName, "a")
        
        while not self.exit and (not self.finished or \
                writtenSentences < self.processedSentences):
            
            thread = self.threads.get(True)
            thread.join()
            result = thread.getResult()
            server = thread.getServer()
            self.servers.put(server, True) 
            tuples = method(result)
            if result == None:
                continue
            if result == False or tuples == False:
                self.exit = True
                print "Last processed line: %d." % self.processedLines
            else:
                outputFile.write(tuples)
            writtenSentences += 1

        outputFile.close()
    
    def signalHandler(self, signum, frame):
        print "Last processed line: %d." % self.processedLines
        self.exit = True

    def run(self, inputFileName, outputFileName, method, startingLine = 0):
        self.cons_thread = threading.Thread(target=self.consumer, 
                args=(outputFileName, startingLine, method))
        self.cons_thread.start()

        signal.signal(signal.SIGINT, self.signalHandler)
        self.producer(inputFileName, startingLine)

        self.cons_thread.join()

if __name__ == "__main__":

    cores = 0
    starting_port = 0
    input_file = ""
    starting_line = 0
    output_file = ""
    host_name = "127.0.0.1"
    method = methods.allrels
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:s:d:n:", ["help", "method",
            "starting_line", "output_dir", "--host_name"])
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
        if o in ("-n", "--host_name"):
            host_name = a
    if len(args) != 3:
        print "Wrong number of parameters."
        print "for help use --help"
        sys.exit(1)
    cores = int(args[0])
    starting_port = int(args[1])
    input_file = args[2]

    if output_file == "": output_file = input_file + ".rel"

    #print "Cores: ", cores
    #print "Staring port: ", starting_port
    #print "Input file: ", input_file
    #print "Starting line: ", starting_line
    #print "Output dir: ", output_dir
    #print "Host name: ", host_name

    generator = Generator(cores, starting_port, host_name)
    generator.run(input_file, output_file, method, starting_line)
    sys.exit(0)

