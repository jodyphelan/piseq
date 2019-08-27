import os
from github import Github

g = Github(os.environ["GH_AUTH_TOKEN"])
repo = g.get_repo("jodyphelan/piseq")
pr = repo.get_pull(int(os.environ["CIRCLE_PR_NUMBER"]))

commit = list(pr.get_commits())[0]
print(commit)
commit.create_comment("hello!!!")
