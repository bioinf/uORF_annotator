import os
import tempfile
import atexit


class TemporaryFileManager:
	def __init__(self, tmp_file_name):
		self.name = tmp_file_name

	def delete_tmp_file(self):
		os.remove(self.name)

	def register_at_exit(self):
		atexit.register(self.delete_tmp_file)

	@classmethod
	def create(cls, suffix, tmp_dir=None):
		tmp_file_name = tempfile.NamedTemporaryFile(dir=tmp_dir, suffix=suffix).name
		return cls(tmp_file_name)
