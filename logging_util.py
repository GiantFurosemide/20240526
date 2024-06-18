"""
module to log messages to a file
"""
import datetime

# input message, print date /time and save to log file
def logging(msg,logfile=__file__+".log"):
    """
    input:
    msg: str, message to log
    logfile: str, path to log file, default is __name__+".log"
    
    output:
    None
    """
    date_time = datetime.datetime.now()
    with open(logfile, 'a') as f:
        f.write(f"[{date_time}] : {msg}\n")
    print(f"[{date_time}] : {msg}")


