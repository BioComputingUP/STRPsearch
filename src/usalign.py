import subprocess

# Wrapper for US-align structure aligner
class Usalign:

    def __init__(self, exe_path="USalign"):
        # Just store executable path
        self.exe_path = exe_path

    def __call__(self, file1, file2, cwd_path='.'):
        # Define US-align command arguments
        # US-align supports both PDB and mmCIF directly
        exe_args = [self.exe_path, file1, file2,'-mol','prot', "-outfmt", "-1"]  # TM-score, RMSD, and alignment
        # Run the alignment and return stdout
        return subprocess.run(
            exe_args,
            cwd=cwd_path,
            capture_output=True,
            encoding='utf-8',
            timeout=120
        ).stdout
