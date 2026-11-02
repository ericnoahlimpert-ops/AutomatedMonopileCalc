"""Progress log of a calculation.

Every message is written with a time stamp and flushed at once, so
run_abaqus.bat / run_abaqus_gui.bat (via run_abaqus.ps1) can show the
progress while Abaqus/CAE is running; the print output of Abaqus/CAE is
not shown live in the console.

The log is written to RESULTS_DIR/progress.log and, if set, to the file
in the environment variable MONOPILE_PROGRESS_LOG (set by run_abaqus.ps1).
"""

import os
import time

ENV_VARIABLE = "MONOPILE_PROGRESS_LOG"


def duration(seconds):
    """h:mm:ss"""
    seconds = int(round(seconds))
    return "%d:%02d:%02d" % (seconds // 3600, seconds // 60 % 60,
                             seconds % 60)


class ProgressLog(object):
    def __init__(self, path):
        self.paths = [path] + ([os.environ[ENV_VARIABLE]]
                               if os.environ.get(ENV_VARIABLE) else [])
        for p in self.paths:
            folder = os.path.dirname(p)
            if folder and not os.path.isdir(folder):
                os.makedirs(folder)
            open(p, "w").close()
        self.start = time.time()

    def write(self, message):
        line = "%s  %s" % (time.strftime("%H:%M:%S"), message)
        print(line)
        for p in self.paths:
            try:
                with open(p, "a") as f:
                    f.write(line + "\n")
            except (IOError, OSError):
                pass                        # e.g. folder deleted meanwhile

    def elapsed(self):
        return duration(time.time() - self.start)
