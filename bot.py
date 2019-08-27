import os
from github import Github

g = Github(os.environ["GH_AUTH_TOKEN"])
repo = g.get_repo("jodyphelan/piseq")

pr = repo.get_pull(int(os.environ["CIRCLE_PULL_REQUEST"].split("/")[-1]))

commit = repo.get_commit(sha=os.environ["CIRCLE_SHA1"])

commit.create_comment("hello %s" % os.environ["CIRCLE_SHA1"])
