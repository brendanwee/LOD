import os, subprocess

# This function was taken from a github post and is not my own.
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def InstallPyDev():
    p = subprocess.Popen(["apt-get", "install", "python-dev"], stdout=subprocess.PIPE, shell=True)
    output,err = p.communicate()
    del err
    p.wait()
    if "not found" in output:
        p = subprocess.Popen(["yum", "install", "python-devl"], stdout=subprocess.PIPE, shell=True)
        output,err = p.communicate()
        del err
        p.wait()
    "https://stackoverflow.com/questions/21530588/fatal-error-python-h-no-such-file-or-directory"
    if not which("pip"):
        print "pip not found, installing pip"
        os.popen("curl -O https://bootstrap.pypa.io/get-pip.py")
        os.popen("sudo python get-pip.py --user")
        os.popen("ln -s ~/.local/bin/pip /usr/sbin/")
    if not which("RiboCode"):
        print "RiboCode not found, installing RiboCode"
        _ = os.popen("pip install RiboCode", )


InstallPyDev()