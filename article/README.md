### How to build the PDF ?

We provide a docker file, which saves you the need to install additional packages to build the pdf

In a console (in the article folder), type:

```
sh prepare-docker.sh
./build-with-docker.sh
```

If you get a permission error on Linux, try

```
sudo sh prepare-docker.sh
sudo ./build-with-docker.sh
```

