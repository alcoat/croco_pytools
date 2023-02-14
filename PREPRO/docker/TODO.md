# TODOS

## Docker

Make a smaller docker image :

- Use multi stage docker file => first to compile and other to give docker only with python ?
- see if compiler is not install two times one with apt and other with conda because of wx

Build croco tools for arm64
=> wx tools seems hard to install on this plateform

Remove again root user for better security.

Use conda-lock file to better reproductibility.

## Python tools

```bash
# Format code with
black

# Lint code with
flake8

# Clean import
autoflake --remove-all-unused-imports -i -r .

# Sort import
isort -rc -sl .
# or isort -rc -m 3 .
```

Use autohooks to do that task on commit.
