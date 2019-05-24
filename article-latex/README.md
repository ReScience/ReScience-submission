### [ReScience C](https://rescience.github.io/) article template

This repository contains the Latex (optional) template for writing a ReScience
C article and the (mandatory) YAML metadata file. For the actual article,
you're free to use any software you like as long as you enforce the proposed
PDF style. A tool is available for the latex template that produces latex
definitions from the metadata file. If you use another software, make sure that
metadata and PDF are always synced.

You can also use overleaf with the [provided template](https://www.overleaf.com/read/kfrwdmygjyqw) but in this case, you'll have to enter `metadata.tex` manually.

#### Usage

For a submission, fill in information in
[metadata.yaml](./metadata.yaml), modify [content.tex](content.tex)
and type:

```bash
$ make 
```

This will produce an `article.pdf` using xelatex and provided font.


After acceptance, you'll need to complete [metadata.yaml](./metadata.yaml) with information provided by the editor and type again:

```bash
$ make
```

