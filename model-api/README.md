# model-api

## Developer notes
https://www.webfactory.de/blog/use-ssh-key-for-private-repositories-in-github-actions, which links to [this](https://github.com/marketplace/actions/webfactory-ssh-agent) GitHub Action.
```
ssh-keygen -t rsa -b 4096 -C "git@gitlab.com:nps-ima/uplands-ci.git"
```
with filename ~/.ssh/githubdeploy (see key creation instructions, [here](https://dev.to/n3wt0n/github-deploy-keys-40k5)).
