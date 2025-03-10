# MapsFitting

Machinery to create and fit Isotropic and Anisotropic MC Maps.
These maps could be FullSky ones (not considering Satellite's exposition) or based on DAMPE's Accepted Events Distribution.

This software is able to produce and fit isotropic and anisotropic sky maps !
It uses a TRandom3 random generator, ensuring a very long repetition period.
Two models of maps are generated, with low and high satistics; this permits to evaluate all the differences while fitting these maps and to extract the anisotropy parameters for different energy-bin statistics !

The fit permormed is a Template Fit, where data maps are fitted using template functions (these templates have been generated for the isotropic component and for all anisotropic ones).

    Analitical Templates:
    
    ================
    
    The templates are not biuld using a MC approach but just obtaining TF2's value for a such value of l and b.
    
    In this case of FullSky maps there's no need to normlize again these templates, because being taken from PDFs. The exposizion of the detector needs to be considered in case of realistic and absolute sky maps; in these cases the template maps have been normalized considering the statistics of the data maps. 
    
    Just note that both HS and LS templates will be created; the binning of the HS and LS maps can be changed.
    
    ================
    

The software is also able toc onsider DAMPE's expositure and acceptance, in order to simulate realistic maps.

Here how to execute:

    ./MapsFit (OPTIONAL: const char* template_path)  (OPTIONAL: const char* data_path)

    --> template_path and data_path are OPTIONAL parameters; tehy have to be specified if you don't want to rebuild both data and templates. In case you want use just one of them, you could specify it and the software will automatically rebuild the missing one !

    Template and Data maps generation is a very time-consuming process; this way they are created and saved on disk. If you need just a fit (maybe changing some parameters) there's no need to create new maps, but they relative (or absolute) path could be given to he software as a parameter (as shown above).

To undertand the physical process and the statistical data analysis one, two different statistics (low and high) will be used, in order to see the differences. All the results will be automatically saved on disk.

The performed fit is a chi2 one.

The software parameters can be easily modified directly in the "simulation.conf" file.

In the "SubmitJobs" diretcory there are both authomatic and manual jobs submitter softwares. The first one (used as standard procedure) automatically compiles all the dependency softwares suing a Python 3 script and a json config file.

Enjoy using!!
