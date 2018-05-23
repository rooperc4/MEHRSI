Habitat Model Abundance Indices estimated from RACE Bottom Trawl Survey (after Rooper and Martin 2012)
================
Chris Rooper
May 14, 2018

Purpose
-------

The purpose of this document and its associated functions is to allow the user to compute annual indices of abundance from survey data using a habitat model based on Rooper and Martin (2012). It is designed to be a 2-stage model, with the first stage setting the biological limits for a species (such as depth ranges from 85-250 m for Pacific Ocean perch or geographical limits on Northern rockfish distribution). In the second stage a iterative modeling process fits CPUE data to habitat variables chosen by the user. Individual habitat variables are constrained to using one of 3 shapes for the fitting (dome, asymptotic or linear). Model fitting is completed using the optim function in R and parameters and variables are removed sequentially until there is no improvement in model AIC. Full details of the modeling and the manuscript can be found at <https://www.st.nmfs.noaa.gov/spo/FishBull/1101/rooper.pdf> and the preceeding papers where the method was developed (Rooper and Martin 2009 - <http://www.int-res.com/articles/meps2009/379/m379p253.pdf>, Rooper 2008 - <http://fishbull.noaa.gov/1061/rooper.pdf>, and Rooper and Zimmermann 2005 - <http://www.int-res.com/articles/meps2005/290/m290p251.pdf>).

------------------------------------------------------------------------

Setup data
----------

### Define species of interest and habitat variables

This code is designed to be run with output from RACEBASE (or other data source) where rows represent hauls to be used in the calculation of the survey index and columns indicate CPUE for species (zero-filled data), habitat and other variables. An example is the data for juvenile Pacific Ocean perch included in the MEHRSI package.

``` r
data("Juvenile_POP_Data")
pander::pandoc.table( Juvenile_POP_Data[1:6,])
```

    ## 
    ## ---------------------------------------------------------------------
    ##  hauljoin    lat     long    year   tdepth   btemp   bdepth   slope  
    ## ---------- ------- -------- ------ -------- ------- -------- --------
    ##   881077    56.14   -135.2   1996    9.4      4.4     291     1.882  
    ## 
    ##   881078    55.92   -135.4   1996    25.4     4.9     275     6.926  
    ## 
    ##   881079    55.93   -135.4   1996    36.4     4.2     369     8.288  
    ## 
    ##   881080    55.93   -134.9   1996    10.4     5.4     180     0.7038 
    ## 
    ##   881081    55.89   -134.6   1996    14.4     8.1      85       0    
    ## 
    ##   881082    55.66   -134.6   1996    20.4     5.5     206     0.4036 
    ## ---------------------------------------------------------------------
    ## 
    ## Table: Table continues below
    ## 
    ##  
    ## -------------------------------------------------
    ##   inverts    ttemp   shrimp    juvenile_POP_CPUE 
    ## ----------- ------- --------- -------------------
    ##      0       11.8    0.01286           0         
    ## 
    ##   0.6885     11.56    1.955          5.518       
    ## 
    ##    0.18      10.77      0              0         
    ## 
    ##   0.8072     10.95   0.05597         8.795       
    ## 
    ##  0.0004213   9.916      0              0         
    ## 
    ##      0       12.05   0.02139        0.3564       
    ## -------------------------------------------------

In this data you have a column for the haul identifier, the midpoint positions of the tow, a year variable indicating the survey year, the depth of the thermocline, the bottom temperature and bottom depth and the seafloor slope. The catch of structure forming invertebrates (coral and sponge), the termperature at the thermocline and the total catch of shrimp are also included as potential habitat variables. The final column is the CPUE of juvenile Pacific Ocean perch (in kg/ha). For this analysis, POP were considered juveniles if they were &lt; 250 mm in length.

PRESENCE-ABSENCE MODEL
----------------------

Int the first stage of the modeling, the idea is to reduce the number of zeros in the data set by looking at the overall presence of the species across variables that effect the distribution of fish on a large scale. For example, Northern rockfish occupy depths shallower than about 200 m and generally occur in the Western and Central Gulf of Alaska, but are absent in Southeast Alaska and Northern BC. So even at depths shawllower than 200 m in Southeast Alaska, Northern rockfish will not occur. This is a different kind of zero observation than a zero observation in for Northern rockfish in the Western GOA at 150 m, since the species will not occur in the region. The initial presence-absence modeling looks at variables that effect distribution on a very large scale (regionally). So for the Rooper and Martin (2012) analyses, we used depth, temperature and the longitude of the bottom trawl (called position in the analyses) as our large scale variables.

We did this by creating a subset of variables (depth, bottom temperature and longitude) that were used for the presence-absence calculations. We also specified the years to include in the analyses (leaving out the 2001 GOA survey that did not sample in the eastern GOA).

``` r
variables1<-c("bdepth","btemp","long")
years_pa<-c(1999,2003,2005,2007,2009)
```

For these presence-absence variables a cumulative distribution function is estimated based on the CPUE weighted varaiables. Then a cutoff value (in this case the defaults are lower and upper 5%) are determined and habitat conditions falling outside those limits receive zero values (as in the fish species will likely never occur here due to large scale constraints) and habitat conditions falling inside these limits receive values of 1 (meaning the fish might occur within these habitat conditions). There potentially be many zeros left in your "present" data and there will be some catches that occur in your "absent" data.

The presence-absence functions used here are cdf\_limits which takes the habitat variables, the observed CPUE and a number of bootstrap replicates. It then calculates the cumulative distribution function for each variables and estimates the upper and lower 5% (by default) of the distribution, and using the specified number of bootstraps (default=1000) estimates an error around the estimate. the second function presence\_absence takes the data set and computes a 1 or 0 for presence or absence based on the cdf\_limits. The imputs to the presence\_absence function are the habitat variables and the limits. It is important to note at this point that the cutoff values (upper and lower 5%) can be changed by editing the function.

``` r
data_pa<-cbind(Juvenile_POP_Data[,(variables1)],Juvenile_POP_Data["juvenile_POP_CPUE"],Juvenile_POP_Data["year"])
data_pa <- data_pa[data_pa$year %in% years_pa,] #Subset those years in the year set that you want to use
limits1<-cdf_limits(data_pa["juvenile_POP_CPUE"],data_pa[,(variables1)],1000,0.05)
pander::pandoc.table(limits1)
```

    ## 
    ## ------------------------------------------------------------------
    ##    &nbsp;     cdf_lower   cdf_upper   cdf_lower_SD   cdf_upper_SD 
    ## ------------ ----------- ----------- -------------- --------------
    ##  **bdepth**     87.42       241.3        5.959          11.97     
    ## 
    ##  **btemp**      3.824       6.369        0.3003        0.05808    
    ## 
    ##   **long**     -165.1      -133.9        2.955          0.1748    
    ## ------------------------------------------------------------------

``` r
pred_pa1<-presence_absence(Juvenile_POP_Data[,(variables1)],limits1,Juvenile_POP_Data["juvenile_POP_CPUE"])
```

    ## 
    ## ------------------------------------------
    ##         &nbsp;           Present   Absent 
    ## ----------------------- --------- --------
    ##  **Predicted Present**     842      1222  
    ## 
    ##  **Predicted Absent**      334      2077  
    ## ------------------------------------------
    ## 
    ## Table: Presence and Absence

In the case of juvenile POP, the two steps are separate largely because we used a subset of the data (without the 2001 survey) to compute the limits, but we want to know the presence or absence for the entire data set. Also, if another method for determining presence or absence, it can be used in the place of these functions, as long as it produces an observation of presence or absence for the model.

CATCH-PER-UNIT-OF-EFFORT MODEL
------------------------------

The second stage of the modeling is to use the habitat variables to fit the CPUE data to an abundance model. For this step another set of habitat variables can be used, or even the same set.

For the juvenile POP modeling we will use bottom temperature, bottom depth, seafloor slope and invertebrate abundance. To tell the model which ones we will be using, we create a vector of names corresponding to the column names in our data set.

``` r
var1<-c("btemp","bdepth","slope","inverts")
```

The equations of the model come in three forms, a dome shaped polynomial curve with 3 parameters (a + b*X + c*X2), a exponential function with two parameters (a*X*exp(-b*X)), and a linear function (a*X) with one parameter. In these equations X is the habitat variable. Initially, the model needs to know what form the equations should take for each variable (3, 2 or 1 respecitively). To do this, you set up a vector that corresponds to the vector of habitat variables with the form for each. Usually, you will want to start with the three parameter version and then interatively remove parameters and variables.

``` r
forms1<-c(3,3,3,3)
```

For computation of survey indices, the column containing the year data must also be identified. Then an object is constructed to provide the initial values for the parameters of the model. In this case we start with all parameters = 0. There is a set of parameters for each of the equations in the model forms and also a single parameter for each year in the model.

``` r
year1<-Juvenile_POP_Data["year"]
par1<-rep(0,(sum(forms1)+length(unlist(unique(year1)))))
```

Next we set up the CPUE data. In the case of juvenile POP (and most of the other survey species) the data are not normally distributed. We used a log-transformation with thte addition of 1/2 the minimum positive catch to come closest to normality. Other data transformations can be used here. For convenience sake we also attached the new log-CPUE back to the data frame.

``` r
ob_CPUE1<-Juvenile_POP_Data$juvenile_POP_CPUE
minob<-min(subset(ob_CPUE1,ob_CPUE1>0))*.5
ob_CPUE1<-log(ob_CPUE1+minob)
Juvenile_POP_Data["lnjpop"]<-ob_CPUE1
```

Finally, we specify how many non-parametric bootstraps of the data will be used to compute the error around the survey index. In this case (and by default) it is 500.

``` r
boot1<-500
```

At this point you are ready to fit the model. The model is fit using the iteration function. It requires the list of names of habitat variables to be used in the CPUE model (vars1), the initial forms of the model (forms1), the observed CPUE of the species (ob\_CPUE1), the predicted presence or absence from the presence-absence model (pred\_pa1), the starting values for parameters (par1) and the variable indicating the year of the survey (year1). The model is fit using a call to the optimiz function and then iteratively it reduces the number of parameters and the number of variables in the model until the lowest AIC is attained. It returns a number of variables in a list. In this case we will return the predicted CPUE and the residuals and attach them to our data frame.

``` r
jpopmodel<-iteration(Juvenile_POP_Data[,(var1)],forms1,ob_CPUE1,pred_pa1,par1,year1)
jpopdata_out<-data.frame(Juvenile_POP_Data,pred_pa1,best_PCPUE=jpopmodel$best_PCPUE,resids=jpopmodel$resids)
head(jpopdata_out)
```

<img src="Raw_data.png" width="2100" /><img src="VariableRelationships.png" width="2100" /><img src="ObservedPredicted.png" width="2100" /><img src="QQplot.png" width="2100" />

Next we bootstrap some errors by resampling the data to produce the index and associated variability.

``` r
index_error<-boot_survey_error(jpopmodel,boot1,year1,Juvenile_POP_Data[,(var1)],pred_pa1,ob_CPUE1)
```

Finally, we take all the pieces and put them together to calculate the index of abundance. For this function, it is assumed you have used the log-transformation with the addition of the minimum observation. For other transformations, a modification of this function is needed.

``` r
JpopIndex<-index_calc(index_error[[2]],jpopmodel[["best_model"]],jpopmodel[["best_parameters"]],Juvenile_POP_Data[,(var1)],year1,minob,boot1)
```

<img src="IndexPlot.png" width="2100" />

Look for any significant spatial patterns in the residuals, by fitting a kriging model and adding back into the predicted CPUE

``` r
spatial_resids(Juvenile_POP_Data$long,Juvenile_POP_Data$lat,jpopdata_out$resids,jpopdata_out$best_PCPUE,ob_CPUE1,region="GOA")
```

<img src="spatial_resids.png" alt="Spatial patterns and correlations among residuals" width="2100" />
<p class="caption">
Spatial patterns and correlations among residuals
</p>
