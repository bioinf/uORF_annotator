import os
import tempfile
import atexit


class TemporaryFileManager:
	def __init__(self, tmp_file_name):
		self.name = tmp_file_name

	def delete_tmp_file(self):
		if os.path.exists(self.name):
			os.remove(self.name)

	def register_at_exit(self):
		atexit.register(self.delete_tmp_file)

	@classmethod
	def create(cls, suffix, tmp_dir=None, delete_at_exit=True):
		tmp_file = tempfile.NamedTemporaryFile(dir=tmp_dir, suffix=suffix, delete=False)
		tmp_file.close()  # Close the file to ensure it's created on disk
		tmp_file_manager = cls(tmp_file.name)
		if delete_at_exit:
			tmp_file_manager.register_at_exit()
		return tmp_file_manager
