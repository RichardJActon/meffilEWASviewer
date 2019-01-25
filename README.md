# meffil EWAS viewer

A viewer for [meffil](https://github.com/perishky/meffil/wiki) EWAS objects

...detailed description...

## Installation

``` r
devtools::install_github("RichardJActon/meffilEWASviewer")
```

## Run the Shiny app

There's only one exported function in the package and it runs the Shiny app:

``` r
meffilEWASviewer::ewas_viewer()
```

This will start the Shiny application running locally.

---

# Hosting this app on your own shiny server

Shiny Server and RStudio Connect expect either a ui.R and a server.R or an app.R file. We’re running our application out of a package with none of this, so we won’t be able to publish it until we fix this problem.

The solution we used is to create a directory called ‘shinyApp’ inside the inst directory of the package. 

The name ‘shinyApp’ was chosen for consistency with Shiny Server which uses a ‘shinyApps’ directory if a user is allowed to serve applications from their home directory.

Inside this directory there is a single ‘app.R’ file with the following line in it:

```
meffilEWASviewer::ewas_viewer()
```

If your shiny server is configured for Apps served [from a users' home directory](http://docs.rstudio.com/shiny-server/#host-per-user-application-directories) (AKA self-publishing) then you can host your own version of the app as described below. If your shiny server is deployed centrally with apps located in `/srv/shiny-server` (usually) your sysadmin can follow similar step to those below.

Shiny Server (and Pro) expect apps to be found in a directory called ‘ShinyApps’, within the users home directory. This means that if we install a Shiny app in a package the final location of the app directory will be inside the installed package, not in the ShinyApps directory. In order to work around this, you can create a link from where the app is expected to be, to where it actually is within the installed package structure.

So you would do something like this in a terminal session:

```
# make sure we’re in our home directory
cd
# change into the shinyApps directory
cd shinyApps
# create a link from our app directory inside the package
ln -s /home/sellorm/R/x86_64-pc-linux-gnu-library/3.4/meffil_ewas_viewer/shinyApp ./ewas-viewer
```

Note: The path you will find your libraries in will differ from the above. Check by running `.libPaths()[1]` and then `dir(.libPaths()[1])` to see if that’s where your packages are installed.

Once this is done, the app should be available at ‘http://<server-address>:3838/<your_username>/ewas-viewer’ and can be updated by updating the installed version of the package. Update the package and the updates will be published via Shiny Server straight away.

# References  / acknowledgements

[meffil](https://github.com/perishky/meffil/wiki)

Min, Josine, Gibran Hemani, George Davey Smith, Caroline L Relton, and Matthew Suderman. 2017. “Meffil: efficient normalisation and analysis of very large DNA methylation samples.” bioRxiv 44 (0): 125963. doi:[10.1101/125963](https://doi.org/10.1101/125963).


This app was created with the help of this guide:

https://www.r-bloggers.com/packaging-shiny-applications-a-deep-dive/

and begun using this template repo: https://github.com/mangothecat/shinyAppDemo 
(linked to form the above guide)

The how to install this app on your own shiny server was cribbed from this guide.
