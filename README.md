# MapsFitting

Machinery to create and fit Isotropic and Anisotropic MC Maps.
These maps could be FullSky ones (not considering Satellite's exposition) or based on DAMPE's Accepted Events Distribution.

This software is able to produce and fit isotropic and anisotropic sky maps !
It uses a TRandom3 random generator, ensuring a very long repetition period.
Two models of maps are generated, with low and hist satistic; thsi permits to evaluate all the differences while fitting these maps ! Obvoiusly HS Tempalte Fits are better and will be used to obtain phisical informations from the data analysis process.

The fit permormed is a Template Fit, where data maps are fitted using template functions (these templates have been generated for the isotropic component and for all anisotropic ones).

Into the repo, however, two different type of Templates are proposed, MC os Anaitical ones.

    MC Templates:

    ================
    This is absolutely the most time consuming way to proceed; in case of HS templates this process could be really long !
    
    In this case Templates are generated using a TRandom3 engine, using PDF's base functions. To use this templates in the fit functions, they have to be proprly normalized (until an alpha version, the correct normalization procedure should be already found !)
    
    ================


    Analitical Templates:
    
    ================
    
    This is a very smart way to proceed, and is the standard one.
    In this time computation time is really negliable respect to the upper case. In this case, infact, templates are not biuld using a MC approach but just obtaining PDF's (FT2) value for a sucg value of l and b.
    
    In this case there's no need to normlize again these templates, because being taken from PDF0s, they're just ok and can be direcly used.
    
    Just note that both HS and LS templates will be created; in case of same binning for HS and LS maps this is redundant, because the final template maps are the same and they don't depend on the satistics used (at difference of the MC case). So this way to proceed has been conserved because it doesnt't imact on the software performace and permits to easily obtain templates with different binning ! Maybe it could be removed in different versions.
    
    ================
    

Until now this software produces (and fits) just All Sky Maps, without considering for DAMPE's expositure and acceptance. An update will be done soon, considering also these important effects.

Here how to execute:

    ./MapsFit (OPTIONAL: const char* template_path)  (OPTIONAL: const char* data_path)

    --> template_path and data_path are OPTIONAL parameters; tehy have to be specified if you don't want to rebuild both data and templates. In case you want use just one of them, you could specify it and the software will automatically rebuild the missing one !

    Template and Data maps generation is a very time-consuming process; thsi way they are created and saved on disk. If you need just a fit (maybe changing some parameters) there's no need to create new maps, but they relative (or absolute) path could be given to he software as a parameter (as shown above).

To better undertand the physical process and the statistical data analysis one, two different statistics (low and high) will be used, in order to see the differences. All the results will be automatically saved on disk.

The performed fit is a chi2 one.

Into "engine.cpp" file you could choice between normalized data maps (or not, that's the standard method to proceed).

Enjoy using!!
