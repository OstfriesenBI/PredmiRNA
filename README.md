# PredmiRNA
BTW this is a [markdown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#links) (.md) file, it allows for quick and easy rich text formatting, that can easily be read text only and converted to other formats.

## Contributing
[Git tutorial](https://rogerdudler.github.io/git-guide/index.html)
### Installing conda
conda or miniconda with Python3 is [required](https://conda.io/docs/user-guide/install/index.html). [Installer](https://conda.io/miniconda.html)
If you want to use anything from the local conda enviroment, you have to activate the conda enviroment or add it the binaries to the system path:
``` sh
source ~/miniconda3/bin/activate 
``` 
Deactivate it
``` sh
source ~/miniconda3/bin/deactivate
```   
### Adding those to the system path
On the linux cluster of the HS Emden/Leer we don't have acces to the .bashrc file, but it loads .bash_aliases, so a little hack after the installation completes:
This overwrites the visibility of some system versions with the ones installed by the installer. 
``` sh  
echo 'export PATH=~/miniconda3/bin:$PATH' >> ~/.bash_aliases
```
Make sure you have a git client installed or the git binary available. For easy acces [add a ssh key to github](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/), so that the storage of the github password is not needed.

### Installing snakemake
Using conda:
``` sh
conda install -c bioconda -c conda-forge snakemake
```
Otherwise: Global installation with pip/easyinstall
``` sh
easy_install3 snakemake
#or
pip3 install snakemake
```

