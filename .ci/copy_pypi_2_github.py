# stdlib
import json
import os
import pathlib
import sys
import tempfile
import urllib.parse

# 3rd party
import github
import requests


def get_pypi_releases(project_name):

	pypi_releases = {}

	# Parse PyPI data
	r = requests.get(f"https://pypi.org/pypi/{project_name}/json")
	if r.status_code != 200:
		print(f"Unable to get package data from PyPI for '{project_name}'", file=sys.stderr)

	else:
		pkg_info = json.loads(r.content)

		for release, release_data in pkg_info["releases"].items():
			release_urls = []
			for file in release_data:
				release_urls.append(file["url"])
			pypi_releases[release] = release_urls

	return pypi_releases


def update_github_release(repo, tag_name, release_name, release_message):
	current_assets = []

	try:
		release = repo.get_release(tag_name)

		# Update existing release
		release.update_release(name=release_name, message=release_message)

		# Get list of current assets for release
		for asset in release.get_assets():
			current_assets.append(asset.name)

	except github.UnknownObjectException:

		# Create the release
		release = repo.create_git_release(tag=tag_name, name=release_name, message=release_message)

	return release, current_assets


def get_file_from_pypi(url, tmpdir):
	filename = pathlib.PosixPath(urllib.parse.urlparse(url).path).name

	r = requests.get(url)
	if r.status_code != 200:
		print(f"Unable to download '{filename}' from PyPI. Skipping.", file=sys.stderr)
		return False

	(tmpdir / filename).write_bytes(r.content)

	return True


def copy_pypi_2_github(g, repo_name, github_username, *, release_message='', pypi_name=None):
	repo_name = str(repo_name)
	github_username = str(github_username)

	if not pypi_name:
		pypi_name = repo_name
	pypi_name = str(pypi_name)

	pypi_releases = get_pypi_releases(pypi_name)

	repo = g.get_repo(f"{github_username}/{repo_name}")
	print(repo.name)

	with tempfile.TemporaryDirectory() as tmpdir:
		tmpdir = pathlib.Path(tmpdir)

		for tag in repo.get_tags():
			version = tag.name.lstrip("v")
			if version not in pypi_releases:
				sys.stdout.flush()
				print(f"No PyPI release found for tag '{tag.name}'. Skipping.", file=sys.stderr)
				sys.stderr.flush()
				continue

			print(f"Processing release for {version}")
			release_name = f"Version {version}"
			release_message += f"""
Automatically copied from PyPI.
https://pypi.org/project/{pypi_name}/{version}
"""

			release, current_assets = update_github_release(repo, tag.name, release_name, release_message)

			# pprint(pypi_releases[version])
			# Copy the files from PyPI

			for pypi_url in pypi_releases[version]:
				filename = pathlib.PosixPath(urllib.parse.urlparse(pypi_url).path).name
				# print(filename)

				if filename in current_assets:
					sys.stdout.flush()
					print(f"File '{filename}' already exists for release '{tag.name}'. Skipping.", file=sys.stderr)
					sys.stderr.flush()
					continue

				if get_file_from_pypi(pypi_url, tmpdir):
					print(f"Copying {filename} from PyPi to GitHub Releases")
					release.upload_asset(str(tmpdir / filename))
				else:
					continue


if __name__ == '__main__':
	gh_token = os.environ.get("GITHUB_TOKEN")
	if not gh_token:
		sys.stdout.flush()
		print("Please supply a GitHub token via the environment variable `GITHUB_TOKEN`.", file=sys.stderr)
		sys.stderr.flush()
		sys.exit(1)

	g = github.Github(gh_token)

	rate = g.get_rate_limit()
	remaining_requests = rate.core.remaining
	print(rate)

	github_username = "domdfcoding"

	copy_pypi_2_github(g, "PyMassSpec", "domdfcoding", pypi_name="PyMassSpec")

	rate = g.get_rate_limit()
	used_requests = remaining_requests - rate.core.remaining
	print(f"Used {used_requests} requests. {rate.core.remaining} remaining. Resets at {rate.core.reset}")
