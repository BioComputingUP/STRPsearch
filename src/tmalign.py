import subprocess


# Wrapper for TM-align structure aligner
class Tmalign:

    def __init__(self, exe_path="TMalign"):
        # Just store executable path
        self.exe_path = exe_path

    def __call__(self, pdb1, pdb2, cwd_path='.'):
        # Define arguments
        exe_args = [self.exe_path, pdb1, pdb2, "-fast", "-a", "T"]
        # Run process
        return subprocess.run(exe_args, cwd=cwd_path, capture_output=True, encoding='utf-8', timeout=120).stdout