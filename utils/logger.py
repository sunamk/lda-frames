#
#copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

import logging
import threading
import sys

log = logging.getLogger("lda-frames")

def setLogFile(filename):
    global log
    hdlr = logging.FileHandler(filename)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    log.addHandler(hdlr)

def setLevel(level):
    global log
    log.setLevel(level)

def DBG(message):
    global log
    log.debug("Thread[%d]: %s" % (threading.current_thread().ident, str(message)))

def INFO(message):
    global log
    sys.stderr.write(message + "\n")
    log.info("Thread[%d]: %s" % (threading.current_thread().ident, str(message)))

def WARN(message):
    global log
    sys.stderr.write(message + "\n")
    log.warning("Thread[%d]: %s" % (threading.current_thread().ident, str(message)))

def ERR(message):
    global log
    sys.stderr.write(message + "\n")
    log.error("Thread[%d]: %s" % (threading.current_thread().ident, str(message)))


