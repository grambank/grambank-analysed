# grambank-analysed
This is a git repository with R scripts for analysis of Grambank data for the following paper:

```
Skirgård, H., Haynie, H. J., Hammarström, H., Blasi, D. E., Collins, J., Latarche, J., Lesage, J., Weber, T., Witzlack-Makarevich, A., Passmore, S., Maurits, L., Dunn, M., Reesink, G., Singer, R., Bowern, C., Epps, P., Hill, J., Vesakoski, O., Robbeets, M., Abbas, K., Auer, D., Bakker, N., Barbos, G., Borges, R., Danielsen, S., Dorenbusch, L., Dorn, E., Elliott, J., Falcone, G., Fischer, J., Ghanggo Ate, Y., Gibson, H., Göbel, H., Goodall, J., Gruner, V., Harvey, A., Hayes, R., Heer, L., Herrera Miranda, R., Hübler, N., Huntington-Rainey, B., Ivani, J., Johns, M., Just, E., Kashima, E., Kipf, C., Klingenberg, J., König, N., Koti, K., Kowalik, R., Krasnoukhova, O., Lindvall, N., Lorenzen, M., Lutzenberger, H., Martins, T., Mata German, C., Meer, S., Montoya Samamé, J., Müller, M., Muradoglu, S., Neely, K., Nickel, J., Norvik, M., Oluoch, C. A., Peacock, J., Pearey , I., Peck, N., Petit, S., Pieper, S., Poblete, M., Prestipino, D., Raabe, L., Raja, A., Reimringer, J., Rey, S., Rizaew, J., Ruppert, E., Salmon, K., Sammet, J., Schembri, R., Schlabbach, L., Schmidt, F., Skilton, A., Smith, W. D., Sousa, H., Sverredal, K., Valle, D., Vera, J., Voß, J., Witte, T., Wu, H., Yam, S., Ye 葉婧婷, J., Yong, M., Yuditha, T., Zariquiey, R., Forkel, R., Evans, N., Levinson, S. C., Haspelmath, M., Greenhill, S. J., Atkinson, Q. D. and Gray, R. D. (submitted) Grambank data reveal global patterns in the structural diversity of the world's languages.
```
This repository also includes analysis and/or plots which requires data from glottolog, AUTOTYP and WALS. All data, including grambank data, is made available in this repository through the use of git submodules.

You can run almost all the scripts in this repository on your personal computer, with the exception of the BRMS analysis and also possibly the INLA-analysis (depending on if it's possible to install INLA on your machine). You can read more about the possible Makefile rules you can run [here](https://github.com/grambank/grambank-analysed/blob/main/R_grambank/README.md#files).

## Git submodules
This Git repository contains git submodules. That means that this repository is linked to other git repositories in a principled way. In this instance this repository has git submodules for the following repostiroeies: grambank-cldf, AUTOTYP-data, glottolog-cldf and WALS.

If you want to run scripts in this repository on your machine, it is necessary not only to clone this repository but also after cloning to run:

`git submodule update --init`

This command will initialise and update the git submodules appropriately. Note that this includes data from grambank-cldf, so no script will run without initalising the git submodules.

Be aware that we have checked out glottolog-cldf with the tag 4.4 in particular, because this is the version of Glottolog that was used to generate grambank-cldf. Checking out later versions of glottolog will cause discontinuity errors.

You can read more about git submodules [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules#_cloning_submodules).
