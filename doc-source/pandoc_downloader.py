# -*- coding: utf-8 -*-
"""
Script to download pandoc on Linux

Based on https://github.com/bebraw/pypandoc
MIT Licensed
"""

# stdlib
import os
import os.path
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.request import urlopen


def _get_pandoc_urls(version="latest"):
	"""Get the urls of pandoc's binaries
	Uses sys.platform keys, but removes the 2 from linux2
	Adding a new platform means implementing unpacking in "DownloadPandocCommand"
	and adding the URL here

	:param str version: pandoc version.
		Valid values are either a valid pandoc version e.g. "1.19.1", or "latest"
		Default: "latest".

	:return: str pandoc_urls: a dictionary with keys as system platform
		and values as the url pointing to respective binaries

	:return: str version: actual pandoc version. (e.g. "lastest" will be resolved to the actual one)
	"""
	# url to pandoc download page
	url = f"https://github.com/jgm/pandoc/releases/{('tag/' if version != 'latest' else '')}{version}"
	# read the HTML content
	response = urlopen(url)
	content = response.read()
	# regex for the binaries
	regex = re.compile(r"/jgm/pandoc/releases/download/.*\.(?:msi|deb|pkg)")
	# a list of urls to the bainaries
	pandoc_urls_list = regex.findall(content.decode("utf-8"))
	# actual pandoc version
	version = pandoc_urls_list[0].split('/')[5]
	# dict that lookup the platform from binary extension
	ext2platform = {'msi': 'win32', 'deb': 'linux', 'pkg': 'darwin'}
	# parse pandoc_urls from list to dict
	# py26 don't like dict comprehension. Use this one instead when py26 support is dropped
	pandoc_urls = {ext2platform[url_frag[-3:]]: ("https://github.com" + url_frag) for url_frag in pandoc_urls_list}
	# pandoc_urls = dict((ext2platform[url_frag[-3:]], ("https://github.com" + url_frag)) for url_frag in pandoc_urls_list)
	return pandoc_urls, version


def _make_executable(path):
	mode = os.stat(path).st_mode
	mode |= (mode & 0o444) >> 2  # copy R bits to X
	print("* Making %s executeable..." % (path))
	os.chmod(path, mode)


def download_pandoc(url=None, targetfolder=None, version="latest"):
	"""Download and unpack pandoc

	Downloads prebuild binaries for pandoc from `url` and unpacks it into
	`targetfolder`.

	:param str url: URL for the to be downloaded pandoc binary distribution for
		the platform under which this python runs. If no `url` is give, uses
		the latest available release at the time pypandoc was released.

	:param str targetfolder: directory, where the binaries should be installed
		to. If no `targetfolder` is give, uses a platform specific user
		location: `~/bin` on Linux, `~/Applications/pandoc` on Mac OS X, and
		`~\\AppData\\Local\\Pandoc` on Windows.
	"""
	# get pandoc_urls
	pandoc_urls, _ = _get_pandoc_urls(version)

	pf = sys.platform

	# compatibility with py3
	if pf.startswith("linux"):
		pf = "linux"
		if platform.architecture()[0] != "64bit":
			raise RuntimeError("Linux pandoc is only compiled for 64bit.")

	else:
		raise RuntimeError("Only Linux is supported.")

	if url is None:
		url = pandoc_urls[pf]

	filename = url.split("/")[-1]
	if os.path.isfile(filename):
		print(f"* Using already downloaded file {filename}")
	else:
		print(f"* Downloading pandoc from {url} ...")
		# https://stackoverflow.com/questions/30627937/tracebaclk-attributeerroraddinfourl-instance-has-no-attribute-exit
		response = urlopen(url)
		with open(filename, 'wb') as out_file:
			shutil.copyfileobj(response, out_file)

	if targetfolder is None:
		targetfolder = "~/bin"
	targetfolder = pathlib.Path(targetfolder).expanduser()

	# Make sure target folder exists...
	if not targetfolder.is_dir():
		targetfolder.mkdir(parents=True)

	print(f"* Unpacking {filename} to tempfolder...")

	tempfolder = pathlib.Path(tempfile.mkdtemp())
	cur_wd = os.getcwd()
	filename = os.path.abspath(filename)
	try:
		os.chdir(tempfolder)
		cmd = ["ar", "x", filename]
		# if only 3.5 is supported, should be `run(..., check=True)`
		subprocess.check_call(cmd)

		dir_listing = set(os.listdir(tempfolder))
		if "data.tar.gz" in dir_listing:
			cmd = ["tar", "xzf", "data.tar.gz"]
		elif "data.tar.xz" in dir_listing:
			cmd = ["tar", "xJf", "data.tar.xz"]
		elif "data.tar.bz" in dir_listing:
			cmd = ["tar", "xjf", "data.tar.bz"]
		else:
			raise FileNotFoundError(f"'data' archive not found. Files in the download are:\n{dir_listing}")

		subprocess.check_call(cmd)
		src = tempfolder / "usr" / "bin" / "pandoc"
		dst = targetfolder / "pandoc"

		if not dst.parent.is_dir():
			dst.parent.mkdir(parents=True)

		print(f"* Copying 'pandoc' to {str(targetfolder)!r} ...")
		shutil.copyfile(src, dst)
		_make_executable(dst)

		src = tempfolder / "usr" / "share" / "doc" / "pandoc" / "copyright"
		dst = targetfolder / "copyright.pandoc"
		print(f"* Copying copyright to {str(targetfolder)!r} ...")
		shutil.copyfile(src, dst)
	finally:
		os.chdir(cur_wd)
		shutil.rmtree(tempfolder)


if __name__ == '__main__':
	download_pandoc()
