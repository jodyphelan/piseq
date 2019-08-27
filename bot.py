import os
from github import Github

g = Github(os.environ["GH_AUTH_TOKEN"])
repo = g.get_repo("jodyphelan/piseq")
print(os.environ["CIRCLE_PULL_REQUEST"])
pr = repo.get_pull(int(os.environ["CIRCLE_PULL_REQUEST"]))

commit = list(pr.get_commits())[0]
print(commit)
commit.create_comment("hello!!!")
