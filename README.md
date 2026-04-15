<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <meta http-equiv="Content-Style-Type" content="text/css">
  <title></title>
  <meta name="Generator" content="Cocoa HTML Writer">
  <meta name="CocoaVersion" content="2575.2">
</head>
<body>
<p class="p1"><span class="s1"></span><br></p>
<h1 style="margin: 0.0px 0.0px 10.0px 0.0px; font: 34.0px 'Helvetica Neue'; color: #262626; -webkit-text-stroke: #262626"><span class="s1">Logistic Growth Modelling Scripts</span></h1>
<p class="p3"><span class="s1">by: Oraya S</span></p>
<h2 style="margin: 0.0px 0.0px 10.0px 0.0px; font: 30.0px 'Helvetica Neue'; color: #262626; -webkit-text-stroke: #262626"><span class="s1">Goal:</span></h2>
<ul class="ul1">
  <li class="li5"><span class="s2"></span><span class="s1">Apply this growth modelling framework to specific pathogen resistance data</span></li>
</ul>
<h2 style="margin: 0.0px 0.0px 10.0px 0.0px; font: 30.0px 'Helvetica Neue'; color: #262626; -webkit-text-stroke: #262626"><span class="s1">Scope + Notes:</span></h2>
<ul class="ul1">
  <li class="li3"><span class="s2"></span><span class="s1">Growth model parameterised in terms on \(\beta\), \(\gamma\), and P_0 - derived from 2 compartment ODE</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">For this version, we predictions are generated from 2000 to 2022. This can be updated.<span class="Apple-converted-space"> </span></span></li>
  <li class="li3"><span class="s2">(updated) </span><span class="s1">Spatial models now use scaled distances instead of distance in Kms - I found that when testing these models, a lot of countries were getting treated like islands from the large distances, which was fine when fitting and estimating trends for countries with data, but when estimating trends for countries with missing data, this resulted in a flat lines (since the model estimates very small spatial correlations for new countries and the values are set to the population mean). Reformatting the model and scaling the distance solved this.</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">Countries are grouped into regions based on definitions from Our World in Data (OWID)</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">spatial models are fitted using cmdstanR because it was computationally more efficient</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">Covariate choice should be selected and changed based on pathogen.</span></li>
</ul>
<h2 style="margin: 0.0px 0.0px 10.0px 0.0px; font: 30.0px 'Helvetica Neue'; color: #262626; -webkit-text-stroke: #262626"><span class="s1">General Flow of Code:</span></h2>
<ul class="ul1">
  <li class="li5"><span class="s2"></span><span class="s1">Standardise data -&gt; fit models -&gt; generate predictions -&gt; generate ensemble predictions</span></li>
</ul>
<h2 style="margin: 0.0px 0.0px 10.0px 0.0px; font: 30.0px 'Helvetica Neue'; color: #262626; -webkit-text-stroke: #262626"><span class="s1">Files:</span></h2>
<ul class="ul1">
  <li class="li3"><span class="s2"></span><span class="s1">load_functions.R : has functions to standardise the data</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">load_packages.R : loads packages you need to run functions in the scripts</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">generatePredictions_functions.R : Has functions to extract predictions in R after models are fitted</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">1_fitModels.R : Fit the models</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">2_extractPredictions.R : Load models, extract predictions, extract model weights used for stacking</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">3_modelStack.R : Load the models, load the stacking weights, generate ensemble predictions</span></li>
</ul>
<p class="p3"><span class="s1">Main Models (Nov 12, 2025):</span></p>
<ul class="ul1">
  <li class="li3"><span class="s2"></span><span class="s1">betaGammaI0_country.stan : country-level random effects</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">betaGammaI0_country_covar.stan : country-level random effects + covariates</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">betaGammaI0_country_AddSource.stan : country and data source-level random effects</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">spatial_cmdStan.stan : spatial-effects</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">spatialcmdStanSource.stan : spatial and data source-level effects</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">spatialcmdStanSourceCovar.stan : spatial and data source-level effects + covariates</span></li>
</ul>
<h2 style="margin: 0.0px 0.0px 10.0px 0.0px; font: 30.0px 'Helvetica Neue'; color: #262626; -webkit-text-stroke: #262626"><span class="s1">Running the Code (simplified)</span></h2>
<ol class="ol1">
  <li class="li3"><span class="s2"></span><span class="s1">Starting with 1_fitModels.R: load the functions and standardise your data. Make sure to set your path. Save your standardised data file as a csv so you can just call it again without having to run the standardising function. If you see errors when loading the functions, run the code line by line in the function file to see what is causing your errors!</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">Run 2_extractPredictions.R: After your fitted models are saved, it’s time to generate some predictions. Now you will load an extra function file called generatePredictions_functions. If you see error when running your code, run everything line by line to see where the errors are coming from. Bottom of the code also generates some plots.</span></li>
  <li class="li3"><span class="s2"></span><span class="s1">Run 3_modelStack.R : Now that you extracted the model weights in the previous step, we will call the file with those model weights and use that to determine the number of samples is allocated to each model to generate ensemble predictions. Plots are generated for the stacked predictions.</span></li>
</ol>
</body>
</html>
