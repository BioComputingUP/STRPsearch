import subprocess


# Wrapper for Foldseek structure aligner
class Foldseek:

    def __init__(self, exe_path="foldseek"):
        # Store executable path
        self.exe_path = exe_path

    def easy_search(self, input_file, db, output_file, temp_dir, columns, eval="1.000E+01", exhaustive="0",
                    aln_type="2", cwd_path='.'):
        """Search one/multiple queries against a target database"""
        # Define arguments
        exe_args = [self.exe_path, "easy-search", input_file, db, output_file, temp_dir, "--format-output",
                    columns, "-e", eval, "--exhaustive-search", exhaustive, "--alignment-type", aln_type]
        # Run process
        return subprocess.run(exe_args, cwd=cwd_path, capture_output=True, encoding='utf-8', timeout=120).stdout

    def createdb(self, input_file, output_file, cwd_path='.'):
        """Create a Foldseek 3Di database from a set of structures"""
        # Define arguments
        exe_args = [self.exe_path, "createdb", input_file, output_file]
        # Run process
        return subprocess.run(exe_args, cwd=cwd_path, capture_output=True, encoding='utf-8', timeout=120).stdout

    def search(self, query_db, target_db, alignment_db, temp_dir, exhaustive="1", cwd_path='.'):
        """Search two Foldseek 3Di databases against each other and save the results as an alignment database"""
        # Define arguments
        exe_args = [self.exe_path, "search", query_db, target_db, alignment_db, temp_dir,
                    "--exhaustive-search", exhaustive]
        # Run process
        return subprocess.run(exe_args, cwd=cwd_path, capture_output=True, encoding='utf-8').stdout

    def createtsv(self, query_db, target_db, alignment_db, output_path, cwd_path='.'):
        """Extract the data in an alignment database in TSV format"""
        # Define arguments
        exe_args = [self.exe_path, "createtsv", query_db, target_db, alignment_db, output_path]
        # Run process
        return subprocess.run(exe_args, cwd=cwd_path, capture_output=True, encoding='utf-8').stdout

