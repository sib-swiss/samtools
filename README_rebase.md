Rebasing to upstream samtools
=============================

```sh
git checkout -b github
git fetch github
git branch --set-upstream-to=github/develop github
git pull
git diff
vi Makefile 
vi bamtk.c 
git diff
vi bamtk.c 
git commit -a -m "merge from github"
git pull
git status
git branch
git checkout develop
git status
git pull
git rebase github
git status
git log
cd ../htslib/
git pull
make clean
make
cd ../samtools/
make clean
make
git push
git branch -D github
```
