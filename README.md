# grambank-analysed
This is a git repository with R scripts for analysis of Grambank data for the following paper:


Skirgård, H., Haynie, H. J.,  Blasi, D. E., Hammarström, H., Collins, J., Latarche, J., Lesage, J., Weber, T., Witzlack-Makarevich, A., Passmore, Chira, A., Dinnage, R., S., Maurits, L., Dinnage, R., Dunn, M., Reesink, G., Singer, R., Bowern, C., Epps, P., Hill, J., Vesakoski, O., Robbeets, M., Abbas, K., Auer, D., Bakker, N., Barbos, G., Borges, R., Danielsen, S., Dorenbusch, L., Dorn, E., Elliott, J., Falcone, G., Fischer, J., Ghanggo Ate, Y., Gibson, H., Göbel, H., Goodall, J., Gruner, V., Harvey, A., Hayes, R., Heer, L., Herrera Miranda, R., Hübler, N., Huntington-Rainey, B., Ivani, J., Johns, M., Just, E., Kashima, E., Kipf, C., Klingenberg, J., König, N., Koti, K., Kowalik, R., Krasnoukhova, O., Lindvall, N., Lorenzen, M., Lutzenberger, H., Martins, T., Mata German, C., Meer, S., Montoya Samamé, J., Müller, M., Muradoglu, S., Neely, K., Nickel, J., Norvik, M., Oluoch, C. A., Peacock, J., Pearey , I., Peck, N., Petit, S., Pieper, S., Poblete, M., Prestipino, D., Raabe, L., Raja, A., Reimringer, J., Rey, S., Rizaew, J., Ruppert, E., Salmon, K., Sammet, J., Schembri, R., Schlabbach, L., Schmidt, F., Skilton, A., Smith, W. D., Sousa, H., Sverredal, K., Valle, D., Vera, J., Voß, J., Witte, T., Wu, H., Yam, S., Ye 葉婧婷, J., Yong, M., Yuditha, T., Zariquiey, R., Forkel, R., Evans, N., Levinson, S. C., Haspelmath, M., Greenhill, S. J., Atkinson, Q. D. and Gray, R. D.  (In prep) "Grambank reveals the importance of genealogical constraints on linguistic diversity and highlights the impact of language loss"

This repository also includes analysis and/or plots which requires data from glottolog, AUTOTYP and WALS. All data, including grambank data, is made available in this repository through the use of git submodules.

You can run almost all the scripts in this repository on your personal computer, with the exception of the BRMS analysis and also possibly the INLA-analysis (depending on if it's possible to install INLA on your machine). You can read more about the possible Makefile rules you can run [here](https://github.com/grambank/grambank-analysed/blob/main/R_grambank/README.md).

## Git submodules
This Git repository contains git submodules. That means that this repository is linked to other git repositories in a principled way. In this instance this repository has git submodules for the following repostiroeies: grambank-cldf, AUTOTYP-data, glottolog-cldf and WALS.

If you want to run scripts in this repository on your machine, it is necessary not only to clone this repository but also after cloning to run:

`git submodule update --init`

This will initialise the submodules in this repos, i.e. 
* grambank (the Grambank data itself)
* glottolog-cldf
* autotyp-data
* wals

This command will initialise and update the git submodules appropriately. Note that this includes neessary data from grambank-cldf, so no script will run without initalising the git submodules.

The grambank data repos contains within it submodules of its own. One of these, `raw/Grambank`, refers to a private repos. If you do not have access to that repos, you will not be able to clone it as a submodule. For this reason, it is not advisable to run `--recursive` when you're updating the submodules as it'll mean you will fail at this cloning of the private repos.

Be aware that we have checked out glottolog-cldf with the tag `v4.4` in particular, because this is the version of Glottolog that was used to generate grambank-cldf. Checking out later versions of glottolog will cause discontinuity errors.

You can read more about git submodules [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules#_cloning_submodules).

### Before publication
Before the journal article is published, the grambank data is private. This means that the entire repos refering to the grambank data is private, and you will need at least read-access to it to clone it. Furthermore, you will need a access token (similar to a password) to clone this repos to your local machine. You can learn more about this [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token).

### Using grambank-analysed as a submodule in another project
If you are working with grambank data for another project, you might find some of the scripts here useful. For example, there are R-scripts in this repos that takes the grambank data, reduces dialects and binarises it correctly. There are also scripts here that prune the EDGE-tree approrpiately.

If you find that helpful, you can add this repos as a submodule in another project and then call on the scripts herein to do those tasks. If you do that, you will need to initalise not only the grambank-analused submodule but also the submodules of it. This is easiest done by navigating into each dir (`cd grambank-analysed/grambank`) and running `git submodule update --init`.


