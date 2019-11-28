#!/bin/sh

setup_git() {
  git config --global user.email "travis@travis-ci.org"
  git config --global user.name "Travis CI"
}

commit_website_files() {
  git checkout -b gh-pages
  git add --all  
#  git add . *.html
  git commit --message "Travis build: $TRAVIS_BUILD_NUMBER"
}

upload_files() {
  git remote add origin-pages https://${GITHUB_TOKEN}@github.com/${GITHUB_PATH}.git > /dev/null 2>&1
  git push --quiet --set-upstream origin-pages gh-pages --force
}

setup_git
commit_website_files
upload_files
