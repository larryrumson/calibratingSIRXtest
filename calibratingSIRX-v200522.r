#######################################################################################################################
## Calibrating SIRX Poper Source Code
## Author: Michael Halem (c) 2020
## This code is released, as to the user's choice: either under the BSD License, the GPL License, or the LGPL License.
## The intention is to give the code users full freedom to incorporate some or all of the code within their modeling work.
## Note: The R packages incorporated may be under different licenses, so that it is the user's responsibility
## to decide which license can be correctly used when incorporating the code in your project.
##
## If this code should be found useful, the author would appreciate a small citation.
##
## A DOI Link to the paper will be here: [DOI Link]
## The first version of this paper is being released on Medrxiv and is entitled
## "Calibrating an Epidemic Compartment Model to Seroprevalence Survey Data"
## You may communicate to the author any issues or questions on the code at the contact address in the medrxiv paper.
########################################################################################################################


#####################################################################################################
## genFiguresAndTables()
## Generate the charts and tables for the medrxiv paper.
## This also serves as an example of how to run runSIRX3 from the command line.
## This also shows how the runSIRX3 can be run from the command line.
#####################################################################################################

genFiguresAndTables<-function(
                              nItem=NA  ## 1 = generate Figures 1,2,3, and data for Tables 1, 7, 8, and 9
                              )         ## 2 = generate Figure 4 and Table 5.
{
    if (is.na(nItem)) return()

    if (nItem == 1)
        makeAntiBodyCurve() ## Make the antibody curves for figures 1, 2 and 3.

    if (nItem == 2) ##Make Figure 4 and Table 5
    {
        runSIRX3('NYC',nProject=10, nWindow=70,dataSource=nyc,Qprob=0.05, bAuto=T,
                 betaChgFact=0.51,betaAuto=T,tABcal=116,tMax=137,tChg=20,bSolveXDX='XDX', ABtarg=0.199, ABweight=5e3) 
    }
}
   

################################
## A few items the author finds convenient

options(width=160, digits=5)
##C like printf
printf <- function(...) cat(sprintf(...))

################################
## External packages that are needed
## Fetch these before running from CPAN using your R system, i.e. R-Studio

library(deSolve)
library('minpack.lm')


##Constants used 

DateRefPosix = as.POSIXct('2019-12-31 OO:OO', tz='EST')
DateRefDate = as.Date('2019-12-31', tz='EST')
EST='EST'

##Convenience names for function parameters
nyc=2
nycold=3
nys=4
euro=5
szDataSource=c('Unused','NYC','NYC-OLD','NYS','EuroCDC')

#######################################################################################################################
## Read in some antibody tables

nyAB1 = read.csv("data/nyab1.csv")
perfectAB = read.csv("data/perfectAB.csv")
exampleAB=read.csv("data/exampleAB.csv")


########################################################################################################################
##Read in population tables from the data subdirectory
    
getPopTable<-function(
                      datasource=euro ##the Euro CDC data has different naming nomenclatures and no states than JHU data
                      )
{
    worldPop = read.csv('data/worldPop.csv', as.is=T, stringsAsFactors=F)
    usPop = read.csv('data/usPop.csv', as.is=T, stringsAsFactors=F)
    colnames(worldPop)[1] = 'CountryState'
    colnames(usPop)[1] = 'CountryState'

    if (datasource==euro) ##Has spaces
        {
        pop = worldPop ##return country names as is -- EuroCDC data has no states so not needed here
        ##aliases
        rownames(pop) = pop$CountryState
        pop=rbind(pop, c(CountryState = 'US', pop = pop['United States',2]))
        pop=rbind(pop, c(CountryState = 'UK', pop = pop['United Kingdom',2]))             
        pop=rbind(pop, c(CountryState = 'United_Kingdom', pop = pop['United Kingdom',2]))
        pop=rbind(pop, c(CountryState = 'South_Korea', pop = pop['South Korea',2]))
        }
    else ##"Has Dots"
        {
        worldPop$CountryState = gsub(' ','.', worldPop$CountryState)
        worldPop$CountryState = gsub(',','.', worldPop$CountryState)
        usPop$CountryState =    gsub(' ','.',    usPop$CountryState)
        usPop$CountryState =    gsub(',','.',    usPop$CountryState)
        usPop$CountryState = sprintf("US.%s",    usPop$CountryState)
        
        pop = rbind(worldPop, usPop)
        rownames(pop) = pop$CountryState
        pop=rbind(pop, c(CountryState = 'US', pop = pop['United.States',2]))
        pop=rbind(pop, c(CountryState = 'UK', pop = pop['United.Kingdom',2]))
        }

    ##browser()
    vPop = as.numeric(pop[,2])
    names(vPop) = pop[,1]

    ##aliases  

    
    return(vPop)
}


#############################################################################################
## SIRX3 based on simple integration based on past value
## Run a SIRX (version 3) Integration based on input parameters
## Note: This is a deterministic run that provides extrapolated results based on the parameters.
## The least squares type fitting is included later in the paper.
##
## SIRX3 Uses either a simple Euler Integration (for debugging comparison), or a Faster Runge-Kutta rk() method
## Originaly based on SIR-X model after Maier and Brockmann https://doi.org/10.1101/2020.02.18.20024414
## Descretized in integer time steps using simple linear approximation. (a Runge-Kutta method would converge faster)
## The names variables are defined as in Maier et al.
## Note: in Maier et al the alpha is https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology beta.
## Maier's beta is wikipedia's gamma.
## If you use the code I would appreciate a small citation in your work in addition to Ben Maier's (who invented the concept).
## The below initial values are from Maier et al's ex-Hubei. I have a few suggestions on modifications to fit better the
## underlying physics (i.e. reality) of the epidemic.  You may email me at my public email address: consult@rocketpwr.com
## t = Days since start of simulation
## Compartments:
## S = Succeptible
## I = Infected
## R = Recovered (or died)
## X = Cumulative Reported Cases
## Additional outputs are:
## RfromS = Removed directly from susceptible compartment without infection
## RfromI = Removed from Infected Compartment without passing thru Confirmed Cases  
## This version adds to X15 two more inputs: beta0 and tChg
## tChg: time (relative to 0) when public health authorities change social distancing
## beta0: the beta before the modification in social distancing
## beta: the beta after the modification in social distancing
#############################################################################################

sirX3<-function(N  = 10000,  #population
               X0 = 5,        #initial cases
               dX0 = 1,       #initial change in cases per day
               kIX = 0.03125, #rate from I compartment to X compartment
               beta = 0.375,  #Initial rate from S compartment to I compartment on and after the social distancing changeover
               betaChgFact = 1,  #The factor that beta changes AFTER the change over at tChg 
               gamma = 0.125, #rate from I compartment to R compartment
               Tmax=30,       #days to extend the simulation
               tChg = 0,      #day to change over from beta0 to beta -- setting to zero defaults to old sirX15 model
               dt = 0.002,    #integration steps
               useRK4=T,       ##F = use Halem's version of Euler Method, T=use R's Runge-Kutta rk45dp7 method
               bLinear=T,     #if false, use a step function.  If true linear up to tChg
               retX=NA,        #if either retX or retDX set return only dot product of (isX,isDX)*(X,DX)
               retDX=NA,       #NOTE: IF running Antibody test, make extra space in 1st cell for result!!
               tABtest=NA,        ##return the antibody test result to this time in cell #1 
               ABcurve=perfectAB, ##a supplied antibody curve
               errorRetVal=NA  ##a value so that nls doesn't croak if would return a zero or negative
               )
{
    ##printf("running sirX3 beta=%6f  X0=%6f kIX=%6f betaChgFact=%6f\n",
    ##       beta, X0, kIX, betaChgFact)
    
    nObs = round(Tmax/dt)
    nStep = round(1/dt)
    R0 = X0 * gamma/kIX
    I0 = dX0/kIX
    S0 = N - I0 - X0 - R0
    if (is.na(betaChgFact)) betaChgFact=1
    
    ##browser()
        
    t = 0
    dat = data.frame(t, S=S0, I=I0, X=X0, R=R0)
    
    if (!useRK4)
    {
        ##descrete time steps with no correction for integer step sizes (Euler Method)
        for (i in 2:(nObs+1))
        {
            i0 = i - 1
            t = i0 * dt
            betaUse = beta
            if (!bLinear || tChg<=0)
            {if (t>=tChg) betaUse = beta*betaChgFact}
            else
            ##bLinear
            {
                wt1 = min(t/tChg, 1)
                wt0 = 1 - wt1
                betaUse = beta * (wt0 + betaChgFact * wt1)
            }
            S = S0 - betaUse * S0/N * I0 * dt
            I= I0 + (betaUse * S0 * I0/N - gamma * I0 - kIX * I0) * dt
            X = X0 + kIX * I0 * dt
            R = R0 + (gamma * I0) * dt        
            if (i0 %% nStep == 0)
            {
                dat = rbind(dat, c(t, S, I, X, R))
            }
            S0 = S
            I0 = I
            R0 = R
            X0 = X        
        }
    }
    else
    {
        ##Adapted from R example in man rk Runge-Kutta Method
        ##Have to back to back the RK models in first and second half
        ##browser()
        betaUse<-function(t, beta, betaChgFact, tChg, bLinear)
        {
            if (tChg <= 0) return(beta)
            if(!bLinear) { if (t>=tChg) return(beta*betaChgFact) else return(beta) }
            ##else bLinear
            wt1 = min(t/tChg, 1)
            wt0 = 1 - wt1
            return(beta * (wt0 + betaChgFact * wt1))
        }
        sirxModel<-function(t, x, parms) {
            with(as.list(c(parms, x)), {
                dS = -betaUse(t, beta, betaChgFact, tChg, bLinear=bLinear) * S * I/ N
                dI = betaUse(t, beta, betaChgFact, tChg, bLinear=bLinear) * S * I/N - gamma * I - kIX * I
                dX = kIX * I
                dR = gamma * I
                res = c(dS, dI, dX, dR)
                list(res)
            })
        }

        ## The parameters
        parms = c(beta=beta, gamma=gamma, kIX = kIX)
        ## vector of time steps
        if (length(Tmax) > 1)
            times = Tmax
        else
            times = seq(0, Tmax, length=abs(Tmax+1))
        ##times = seq(0, Tmax*nStep, length=nStep*Tmax+1)/nStep
        xstart = c(S=S0, I=I0, X=X0, R=R0)
	names(xstart) = c('S','I','X','R')
	##atol = c(C=0.01, I=0.01, X=0.01, R=0.01)
        dat = rk(xstart, times, sirxModel, parms, method='rk45dp7', maxsteps=1000) ##Dormond-Prince Runge-Kutta adaptive stepsize        
        dat = as.data.frame(dat)
    }

    cumI = N - dat$X - dat$S  ##This is cumulative current and past tense infected, not current infected
    dat = cbind(dat, cumI)

    X = dat$X
    dX = c(NA, X[-1] - X[-length(X)])
    dat = cbind(dat, dX)


    if (!is.na(tABtest))
    { ##Then calculate an antibody test result
        dCumI = c(dat$cumI[1], dat$cumI[-1] - dat$cumI[-dim(dat)[1]])
        i = which(dat$time <= tABtest)
        NtestPop = N - dat$X[tABtest == dat$time]
        if (length(NtestPop) == 0)
        {
            printf("Error in sirX3() NtestPop failure to calc because of bad dat$time\n")
            browser()
        }
        ABtestResult = testPctFromNumInfected(dCumI[i], tABtest - dat$time[i], nPop = NtestPop, curve=ABcurve)
    }

    dat0=dat
    
    if (!is.na(retX) || !is.na(retDX)) ##a vector of 1s and 0s if needed
    {
        if (!is.na(retX) && is.na(retDX)) retDX = rep(0, length(retX))
        if (!is.na(retDX) && is.na(retX)) retX = rep(0, length(retDX))
        if (length(retDX) != length(retX))
        {
            printf("Error in sirX3() length(retX) != length(retDX))\n");
            return(NA);
        }
        
        X = dat$X
        dX = dat$dX
        nX = length(X)
        remainder = length(retX) %% nX

        if (remainder == 0)  
            dat = X*retX + dX*retDX ##This recycles X and dX if shorter than retDX.
        if (remainder == 1)
        {
            ABtestResult = exp(ABtestResult) #Result is exponentiated so when logged pct relative difference is significant
            dat = c(ABtestResult,  X*retX[-1] + dX*retDX[-1]) ##stick the AB test in first place
        }
        if (remainder>1)
        {
            printf("Error in sirX3() length(retDX) not multiple of length(X) [optional +1]\n");
            return(NA);
        }     
        iError = is.na(dat) | dat<=0
        if (length(errorRetVal) == length(dat)) dat[iError] = errorRetVal[iError]
        if (length(errorRetVal) == length(dat)) dat[iError] = 1 ##If it pops an error, need to return erroneous value
        
    }

    ##debugging print...
    ##if (!is.list(dat)) printf("called sirX3(retX[1,end]=%.0f %.0f retDX[1, end]=%.0f %.0f len=%d %d dat[1:3] %f %f %f)\n",
    ##                          retX[1], retDX[1], retX[length(retX)], retDX[length(retDX)], length(retX), length(retDX),  dat[1], dat[2], dat[3])
    ##print(dat)
    ##if (!is.list(dat))
    ##    browser()
    return(dat)
}


#############################################################################################
## Get the nyc data
## The file nyc.csv is included in the paper's data

getNYC<-function()
{    
    dat = read.csv('data/nyc.csv', stringsAsFactors=F)
    colnames(dat) = gsub('_','.', colnames(dat))
    if (any(colnames(dat) == 'DATE.OF.INTEREST'))
        Date = as.Date(dat$DATE.OF.INTEREST, '%m/%d/%y')
    else
        Date = as.Date(rownames(dat), '%m/%d/%y')
    nDat = dim(dat)[1]
    dC = dat$CASE.COUNT            
    dD = dat$DEATH.COUNT
    C = D = 0
    for (i in 1:nDat)
    {
        C[i] = sum(na.rm=T, dC[1:i])
        D[i] = sum(na.rm=T, dD[1:i])
    }
    if (dC[nDat] < 0.1*dC[nDat-1])
    { ## bad last day filter
        dC = dC[-nDat]
        dD = dD[-nDat]
        C = C[-nDat]
        D = D[-nDat]
        Date = Date[-nDat]
        nDat = nDat - 1
    }

    if (dC[nDat] < 0.6 * dC[nDat - 1] && dC[nDat]>0.4 * dC[nDat - 1])
    { ##second of last day 1/2 day adjustment
        dChalf = dC[nDat]
        dC[nDat] = dC[nDat] + dChalf
        C[nDat] = C[nDat] + dChalf
        
    }

    X = rep(NA, nDat)
    R = rep(NA, nDat)

    t = as.numeric(Date - DateRefDate)

    dat = cbind(Date,t,C,D,R,X)
    dat = as.data.frame(dat, stringsAsFactors=F)
    dat$Date = Date
    dat = dat[order(Date),]
    rownames(dat)=1:dim(dat)[1]
    return(dat)
}



#############################################################################################
## Generalized csv fetcher for case time series

getCsvData<-function(fnamePath, bIsCum=T)
{    
    dat = read.csv(fnamePath, stringsAsFactors=F)
    dat$Date = as.Date(dat$Date)
       
    nDat = length(dat$Date)
    
    Date = dat$Date
    
    C = dat$Cases 
    D = rep(NA, nDat)
    R = rep(NA, nDat)
    X = rep(NA, nDat)
    t = as.numeric(Date - DateRefDate)

    if (!bIsCum)
    {
        dC = C
        for (i in 1:length(C))
        {
            C[i] = sum(na.rm=T, dC[1:i])
        }            
            
    }

    dat = cbind(Date,t,C,D,R,X)
    dat = as.data.frame(dat, stringsAsFactors=F)
    dat$Date = Date
    dat = dat[order(Date),]
    rownames(dat)=1:dim(dat)[1]
    return(dat)
}




#############################################################################################
## Read the NYS county based data to return a data frame used for the

getNysData<-function(county)
{   
    county = toupper(county)
    if (county == "NYC")
    {
        dat1 = getNysData("New.York")
        dat2 = getNysData("Queens")        
        dat3 = getNysData("Kings")
        dat4 = getNysData("Bronx")
        dat5 = getNysData("Richmond")

        tt = c(dat1$t, dat2$t, dat3$t, dat4$t, dat5$t)
        tt = sort(unique(tt))

        dat = dat1[1,] ##clone new data frame
        dat = dat[-1,] ##empty
        for (t in tt)
        {
            Date = dat1$Date[dat1$t == t]
            C = sum(dat1$C[dat1$t == t],
                    dat2$C[dat2$t == t],
                    dat3$C[dat3$t == t],
                    dat4$C[dat4$t == t],
                    dat5$C[dat5$t == t])
            nTst = sum(dat1$nTst[dat1$t == t],
                    dat2$nTst[dat2$t == t],
                    dat3$nTst[dat3$t == t],
                    dat4$nTst[dat4$t == t],
                    dat5$nTst[dat5$t == t])
            newRow = as.data.frame(cbind(Date, t, C, nTst))
            newRow$Date = Date
            dat = rbind(dat, newRow)
            dat$Date = as.Date(dat$Date)
        }
        return(dat)        
    }
    
    county = gsub(' ','.',x=county)
    
    dat = read.csv('data/countyNY.csv', stringsAsFactors=F)
    colnames(dat) = gsub(' ', '.', x=colnames(dat))
    ##Date = as.Date(dat$Test.Date)
    dat$County = toupper(gsub(' ','.',x=dat$County))
       
    counties = sort(unique(dat$County))
    ##browser()    
              
    iCounty = which(county == dat$County)
    
    nDat = length(iCounty)
    Date = as.Date(dat$Test.Date[iCounty], format="%m/%d/%Y")
    C = dat$Cumulative.Number.of.Positives[iCounty] 
    ##D = dat$death[iCountry]
    nTst = dat$Cumulative.Number.of.Tests.Performed[iCounty]
    ##R = dat$recovered[iCountry]
    ##X = rowSums(cbind(D,R), na.rm=T)
    t = as.numeric(Date - DateRefDate)

    ##browser()
    dat = cbind(Date,t,C,nTst)
    dat = as.data.frame(dat, stringsAsFactors=F)
    dat$Date = Date
    dat = dat[order(Date),]
    dat = unique(dat)
    rownames(dat)=1:dim(dat)[1]
    return(dat)
}




#############################################################################################
##Pick the better of the multple models
## In the runSIRX3 curve fitting on bootstrap, a number of initial non-linear closed form models
## are initially fitting, approaching the fit from, for example, r=positive growth or
## r= negative growth.  (modelNls1a and modelNls1b respectively.)
## It was found that sometimes approaching from the wrong side
## would fail to converge, so both sides are run.
## Similarly, the model can be formulated in different ways depending on where the
## constants multiplied through, as seen in modelNls1c and modelNls1d 
## After the 4 modelNls1[a,b,c,d] are run, this function picks the one with the best fit
## and then returns the boot strap parameters, r0, Io, Xo, and dXo
## within the structure so that the runSIRX2 can use this to seed the
## SIRX least square fit.
#############################################################################################

pickBetterModel<-function(kIX=NA, model1, model2, model3, model4)
{
    modelNls = NA
    model = list(model1, model2, model3, model4)

    good = rep(NA, 4)
    sigma = rep(NA, 4)
    for (i in 1:4)
        {
        good[i] = !(class(model[[i]]) == 'try-error')
        if (good[i]) sigma[i] = summary(model[[i]])$sigma
        }

    if (!any(good)) ##then all bad, return first bad one
        return(model1)


    minSigma = min(na.rm=T, sigma)

    iMin = min(which(sigma == minSigma))

    modelNls = model[[iMin]]

    coefs = modelNls$m$getPars()

    Xo = NA
    r0 = NA
    Io = NA

    if (names(coefs)[3] == 'A')
    { ##older model: log(C)~ logOrNA(A * sign(r) * exp(r * (t - tMin))+Co )
        r0 = coefs[1]
        Io = coefs[3] * abs(coefs[1])/kIX 
        Xo = coefs[2] + sign(r0)*coefs[3]
        dXo = abs(coefs[3] * r0) 
    }
    
    if (names(coefs)[3] == 'Apct')
    { ##newer pct model: log(C)~ logOrNA(Co * (1 + Apct * sign(r) * exp(r * (t - tMin))))
        r0 = coefs[1]
        Xo = coefs[2] * (1 + sign(r0)*coefs[3])
        Io = Xo/kIX
        dXo = abs(coefs[2] * coefs[3] * coefs[1])
    }
    
    modPars = list(r0=unname(r0), Xo=unname(Xo), Io=unname(Io), dXo = unname(dXo))
    modelNls$modPars = modPars
    
    return(modelNls)
}

#############################################################################################
## logOrNA() return the log unless R's log() returned a NaN, in which case it returns an NA.
## This is because NAs are handled better by the regression algorithms, while NaNs terminate
## with an error code
#############################################################################################

logOrNA<-function(x) {
    y=rep(NA, length(x))
    i=x>=0
    y[i] = log(x[i])
    return(y)
}


#############################################################################################
## runSIRX3 -- This contains the algorithms that fit least squares to calibrate against the
## seroprevalence data.
## Note that 3 indicates this is roughly the 3rd version of this model (for the medrxiv paper).
## The author will explain in greater detail any thing within.
## If there is was a need, this code could be modularized and made into an R package.
## It could also the basis of a small Github project. 
## If you have further interest or questions, please contact michael.halem@becare.net
#############################################################################################


runSIRX3<-function(country, nWindow=7, tMin=NA, tMax=NA, nProject=2, Cnew=NA,
                  dataSource=euro,
                  N=NA,
                  gamma=0.125, kIX = NA,
                  Qprob=0.2, ##probability of detection and quarantine to X compartment before recovery to R compartment
                  bAuto=F,
                  tChg=NA, ##change over time (after Tmin) from beta0 to beta (final). tChg=0 is equivalent to single beta
                  betaChgFact=1,
                  betaAuto=F,
                  bPlotExtra=F,
                  bPrintOutOfContext=F, ##Print extra internal debugging information -- READ THE CODE!!! for documentation
                  bTrace=F,
                  bRetxdat=F,
                  scale=1,
                  bYtotInterval=NULL,
                  tABcal=NA,
                  ABcurve=nyAB1,
                  bSolveXDX = 'X',  #X = cases, DX=change in cases, XDX=both cases and change in cases
                  ABweight=100,     ##The weight, relative to 1, of the AntiBody log target miss: NOT LOG
                  ABtarg=NA
                  )
{
    if (bPrintOutOfContext)
        printf("WARNING: PRINTING OUT OF CONTEXT EXTRA DEBUGGING INFORMATION ON -- READ THE CODE TO UNDERSTAND!!!")
    if (is.na(tChg)) tChg = nWindow
    ##betaChgFact0 = betaChgFact
    ##rm(betaChgFact)
    
    ##Note: Qprob = kIX/(kIX + gamma) (From Ben Meier's paper)
    ##thus: Qk + Qg = k => Qk - k = -Qg => k(Q - 1) = -Qg
    ##thus: k = Qg/(1 - Q)
    Reff = NA; RinfPct=NA; Xinf=NA  ##variables to put NA if can't calculate so printing won't crash
    if (is.na(kIX))
        kIX = Qprob * gamma/(1 - Qprob)
    Qprob = kIX/(kIX + gamma) ##a doublecheck
    
    ##six=makeplotSix(country, nWindow=0, nProj=0, dTexit=7);       
    if (dataSource==nyc) xdata=getNYC()
    if (dataSource==nys) xdata=getNysData(country)
    if (dataSource==nycold) xdata=getCsvData('data/nycold.csv', bIsCum=T)
    xdata$C=scale*xdata$C    

    if (is.na(N))
        { ##fetch N from population lookup if available
            pop = getPopTable(dataSource)
            N = pop[country]
            if (dataSource==nys) N=pop[sprintf("US.NY.%s", country)]
            if (length(N)  == 0) N = NA
            if (length(N)  >  1) N = N[1]
        }
        
    xdata = xdata[(!is.na(xdata$C)),]


    if (!is.na(Cnew)) #if have today's newest data
        {
            ##browser()
            newData = xdata[dim(xdata)[1],]
            newData$Date = newData$Date + 1
            newData$t = newData$t+1
            newData$C = Cnew
            newData[,4:(dim(newData)[2])]=NA
            xdata = rbind(xdata,newData)
        }

    
    nData = dim(xdata)[1]
    six=xdata
    ##six = cbind(six, dD = c(NA, six$D[-1] - six$D[-nData]))
    if (any(colnames(six) == 'nTst')) six = cbind(six, dTst = c(NA, six$nTst[-1] - six$nTst[-nData]))
    six = cbind(six, dC = c(NA, six$C[-1] - six$C[-nData]))
    if (any(colnames(six) == 'nTst')) six = cbind(six, yTot = six$dC/six$dTst)

    if (!is.null(bYtotInterval))
    {
        interval = match(bYtotInterval, six$t)
        six = cbind(six, Cold = six$C)
        six = cbind(six, dCold = six$dC)
        modYtot = lm(log(yTot)~t, six[interval,])
        meanTst = mean(six$dTst[interval])
        six$dC[interval] = exp(predict(modYtot))*meanTst
        for (i in interval) six$C[i] = six$C[i-1] + six$dC[i]         
    }
    
    ##browser()
    
    Cmax = max(na.rm=T, six$C)
    

    if (is.na(tMax)) tMax = max(na.rm=T, six$t)
    if (is.na(tMin)) tMin = tMin = max(tMax - nWindow, min(na.rm=T, six$t))
    
    iData = six$t>=tMin & six$t<=tMax

    ##double lm method -- first find r, then find A and C
    nDC = sum(na.rm=T, iData & six$dC>0)

    if (nDC <= 1)
    {
        printf("Warning: Flatline Cases -- nothing to predict\n")
        rDC = 0; Aref = 0; Tref = tMin
        Cdc = mean(na.rm=T, six$C[iData])
        modelDC = NULL
        modelDC3 = NULL
    }
    if (nDC > 1)
        {
            
            modelDC = lm(log(dC)~t, data=six[iData & six$dC>0,])
            rDC = modelDC$coefficients[2]
            names(modelDC$coefficients)  = c('Intercept      ', 'r  from dC      ')
            modelDC3 = lm(C ~ I(sign(rDC)*exp(rDC*(t-tMin))), data=six[iData,])
            ##modelDC4 = lm((C) ~ I(sign(rDC)*exp(rDC*(t-tMin))), data=six[iData,])
            Cdc = modelDC3$coefficients[1]
            names(modelDC3$coefficients) = c('Co from dC     ', 'A  from dC      ')
            ##Aref = abs(modelDC3$coefficients[2]/modelDC3$coefficients[1])
            Aref = modelDC3$coefficients[2]
            Tref = -log(Aref)/rDC+tMin
            ##Aref = exp(rDC * tMin)
        }

    mySummary<-function(model)
    {
        if (is.null(model) || class(model) == 'try-error')
        {
            if (is.null(model))
                printf("No model\n")
            else
                {
                    printf("Note: nls model not used due to:\t")
                    print(attributes(model)$condition)
                    printf("\n")
                    printf("\tInitial Parameter Estimates: r=%g Cdc=%g Aref=%g\n", rDC, Cdc, Aref)
                    printf("Possible Fix: re-read source-code i.e. source('makeplot.r')\n")
                }
            }
        else
        {
            szOutput = capture.output(print(width=100, digits=5, summary(model)))
            iBlank = (nchar(szOutput) == 0)
            szOutput = szOutput[!iBlank]
            iResidual = szOutput == "Residuals:" | szOutput == "Call:" | szOutput ==  "Coefficients:" | szOutput == "Parameters:"
            szOutput = szOutput[!iResidual]
            iClutter = grepl('Signif. codes', szOutput)
            szOutput = szOutput[!iClutter]             
            printf("%s\n", szOutput[c(1:length(szOutput))])
        }
    }
        

    ##browser()
    ##Aref = min(abs(Aref/Cdc)) ##New for May 12
    if (bPrintOutOfContext) printf("Aref=%g  Cdc=%g  rDC=%g To=%f\n", Aref, Cdc, rDC, Tref)

  
    Cdc = abs(Cdc)  ##consider when flipping sign on Cdc also fipping sign on A and r
    printf("Running 4 nlsLM boostraps...\n")
    if (bTrace) printf("Running nlsLM #1a\n")
    modelNls1a =  try(silent=T,
                      nlsLM(log(C)~ logOrNA(Co * (1 + Apct * sign(r) * exp(r * (t - tMin)))), data=six[iData,],
##                        start=list(r=rDC, Co=-sign(rDC)*(Cdc), A=abs(Cdc)),
                            start=list(r=unname(rDC), Co=unname(Cdc), Apct=0.5),
                            lower=c(r=-10, Co=-1e10, Apct=-1),
                            upper=c(r=10, Co=1e10, Apct=.9999999),               
                            algorithm='LM', trace=bTrace, control=c(warnOnly=T, maxiter=1000, maxfev=1000)))

    ##Try from the other side
    if (bTrace) printf("Running nlsLM #1b\n")

    modelNls1b =  try(silent=T,
                      nlsLM(log(C)~ logOrNA(Co * (1 + Apct * sign(r) * exp(r * (t - tMin)))), data=six[iData,],
                            ##                        start=list(r=rDC, Co=-sign(rDC)*(Cdc), A=abs(Cdc)),
                            start=list(r=unname(-rDC), Co=unname(Cdc), Apct=0.5),
                            lower=c(r=-10, Co=-1e10, Apct=-1),
                            upper=c(r=10, Co=1e10, Apct=.9999999),               
                            algorithm='LM', trace=bTrace, control=c(warnOnly=T, maxiter=1000, maxfev=1000)))

    if (bTrace) printf("Finished Running nlsLM #1a/b\n")


    if (bTrace) printf("Running nlsLM #1c\n")
    modelNls1c =  try(silent=T,
                      nlsLM(log(C)~ logOrNA(A * sign(r) * exp(r * (t - tMin))+Co ), data=six[iData,],
                              ##                        start=list(r=rDC, Co=-sign(rDC)*(Cdc), A=abs(Cdc)),
                            start=list(r=unname(rDC), Co=unname(Cdc), A=unname(Aref)),
                            lower=c(r=-10, Co=-1e10, A=0),
                            upper=c(r=10, Co=1e10, A=1e10),               
                            algorithm='LM', trace=bTrace, control=c(warnOnly=T, maxiter=1000, maxfev=1000)))

        ## changing sign to approach from other side
    if (bTrace) printf("Running nlsLM #1d\n")
    modelNls1d =  try(silent=T,
                      nlsLM(log(C)~ logOrNA(A * sign(r) * exp(r * (t - tMin))+Co ), data=six[iData,],
                            start=list(r=-unname(rDC), Co=unname(Cdc), A=unname(Aref)),
                            lower=c(r=-10, Co=-1e10, A=0),
                            upper=c(r=10, Co=1e10, A=1e10),               
                            algorithm='LM', trace=bTrace, control=c(warnOnly=T, maxiter=1000, maxfev=1000)))     
        if (bTrace) printf("Finished nlsLM #1c/d\n")

    ##printf("...Finished Running 4 nlsLM boostraps, Up to 3 Failures are OK as long as 1 is Good\n")

    
 
    modelNls = pickBetterModel(kIX, modelNls1a, modelNls1b, modelNls1c, modelNls1d)

    nlsErr = class(modelNls) == 'try-error'   
   
    xdat = six
    if (any(colnames(xdat) == 'X')) xdat=xdat[,-which(colnames(xdat) == 'X')]

    ##Add rows for forward projection
    for (i in 1:nProject)
        {
            newData = xdat[dim(xdat)[1],]
            newData$Date = newData$Date + 1
            newData$t = newData$t+1
            newData[,3:(dim(xdat)[2])]=NA
            xdat = rbind(xdat,newData)
        }
    rownames(xdat)=1:(dim(xdat)[1])
    
    
    nSim = dim(xdat)[1]
    if (!is.null(modelDC3))
        xdat=cbind(xdat, Csim= predict(modelDC3, newdata=data.frame(t=xdat$t)))
    else
        xdat = cbind(xdat, Csim=rep(NA, nSim))
    xdat=cbind(xdat, dCsim=c(NA, (xdat$Csim[-1] - xdat$Csim[-nSim]) ))

    ##FOR comparison -- full nls model
    ##Note: the columns for the full nls model are named Cnlm and dCnlm
    if (nlsErr)
        xdat = cbind(xdat, Cnlm = rep(NA, length(xdat$t)))
    else
        xdat=cbind(xdat, Cnlm= exp(predict(modelNls, newdata=data.frame(t=xdat$t))))
    xdat=cbind(xdat, dCnlm=c(NA, (xdat$Cnlm[-1] - xdat$Cnlm[-nSim])))
    rownames(xdat) = 1:(dim(xdat)[1])

    havFwd = F
    modelSIRXerr=T
    modelSIRX=NULL
    modelSIRXauto = NULL
    kIXauto=NA
    modPars=modelNls$modPars
    if (!is.na(N) && !nlsErr)
    {
        ##coefs = modelNls$m$getPars()
        gammaEff = gamma + kIX
        Reff = 1 + modPars$r0/gammaEff

        if (T) ##simplistic reduction of r
        {
            SIRinf<-function(x) {abs(x - (1-exp(-Reff * x)))}
            root = optimize(interval=c(0,1), f=SIRinf, tol=1e-6)
            ##method='Brent', lower=c(x=0), upper=c(x=1), control=c(trace=9))
            RinfPct = NA
            if (root$objective<1e-6)
                RinfPct = root$minimum
            Rinf = N * RinfPct
            Xinf = Rinf * kIX/(gamma + kIX)
        }

        r0 = modPars$r0
        Xo = modPars$Xo
        Io = modPars$Io
        dXo = modPars$dXo
        
        ##browser()
       
        ##xdat = cbind(xdat, Cfwd, dCfwd)        
        if (r0>0) havFwd = 1

        iMin = which(xdat$t == tMin)
        if (any(iMin))
        {
            R0 = Xo * gamma/kIX  #This is recovered R[t=0] is X times the ratio of the drain factors (gamma/kIX)
            S0 = N - Io - Xo - R0 ##this is susceptibles S[t=0]
            ##beta00 = (r0 + gamma + kIX)*N/S0
            beta00 = (r0 + gamma + kIX)

            ##Ro = 1 + r0/gamma
            ##beta00 = gamma * Ro
            ##Calculate an xdat matrix using the simplistic model above
            ###Now calculated in modPars dXo = r0 * Xo
            xdat15nlm = sirX3(N=N,
                            X0 = Xo,
                            ##dX0 = abs(r0) * coefs[3],
                            dX0 = dXo,
                            kIX = kIX,
                            gamma = gamma,
                            beta = beta00,
                            Tmax = nProject + (tMax - tMin),                            
                            )
            if (bPrintOutOfContext) printf("R0=%f Io=%f r0=%f Xo=%f beta00=%f\n", R0, Io, r0, Xo, beta00)
            ##if (Xo<0) Xo=1
            ##if (Xo<0) = Xo=0.01
            
            dTmax = tMax - tMin
            dTABcal = tABcal - tMin

            dXmax = 0.5*N
            dXmin = 0
            ##NB: dX0 is used to calculate Io. 
             
            start=c(X0=Xo, dX0 = dXo, beta=beta00, betaChgFact=betaChgFact)
            
            lower=c(X0=0, dX0=dXmin, beta=0, betaChgFact=0.0)
            upper=c(X0=+N, dX0=dXmax, beta=Inf, betaChgFact=1)               
            names(start) = names(lower) = names(upper) = c('X0', 'dX0', 'beta', 'betaChgFact')
            
            
            if (tChg <= 0)
            {
                tChg = 0
                betaChgFact = 1
            }

            if (!betaAuto)
                {
                    ##betaChgFact = betaChgFact
                    start = start[1:3]
                    lower = lower[1:3]
                    upper = upper[1:3]                    
                }
            ##browser()
            printf("running modelSIRX\n")
            ##try

            six2 = six[iData,]
            six2 = cbind(six2, isC=rep(0, dim(six2)[1]) , isDC=rep(0, dim(six2)[1]))
            if (bSolveXDX == 'X') six2$isC = 1
            if (bSolveXDX == 'DX') six2$isDC = 1
            if (bSolveXDX == 'XDX')
                {
                    six2C = six2DC = six2
                    six2C$isC = 1
                    six2DC$isDC = 1
                    six2 = rbind(six2C, six2DC)
                }
            iDataNoNeg = (six2$isC & six2$C>0) | (six2$isDC & six2$dC>0)
            six2 = six2[iDataNoNeg,]

            if (!is.na(ABtarg) && !is.na(tABcal))
            {
                six2 = rbind(rep(NA, dim(six2)[2]), six2)
                six2$C[1] = exp(ABtarg)
                six2$isC[1]=1
                six2$isDC[1]=0
                six2$Date[1]=as.Date('1900-01-01') ##Because I have na.omit on regression, must fill in stub data
                six2$t[1]=-99
                six2$nTst[1]= six2$dTst[1]=0
                six2$dC[1]=0
                six2$yTot[1]=0
                ##printf('here!!\n')
                ##browser()
                             
            }

            ##browser()
            if (bTrace) printf("Running nlsLM #3\n")
            modelSIRX =  try(nlsLM(log(C*isC + dC*isDC) ~ log(sirX3(  
                                       retX=isC,
                                       retDX=isDC,
                                       tABtest=dTABcal,
                                       ABcurve=ABcurve,
                                       Tmax=dTmax,
                                       kIX = kIX,
                                       N=N,
                                       errorRetVal=(isC*C + isDC*dC),
                                       tChg = tChg, ##Note: couldn't get it to autosolve for tChg
                                       gamma=gamma,  ##gamma and above are known
                                       ##betaChgFact = betaChgFact0,
                                       ##beta0=beta00, tChg=0,  
                                       X0=X0, dX0=dX0, beta=beta, betaChgFact=betaChgFact ##this is being solved for
                                   )),
                                   data=six2,
                                   start=start,
                                   lower=lower,
                                   upper=upper,
                                   weights = c(ABweight, rep(1, dim(six2)[1]-1)),
                                   algorithm='LM', trace=bTrace, control=list(warnOnly=T, maxiter=1000, maxfev=1000),
                                   na.action = na.omit
                                   ))
            ##browser()
            if (bTrace) printf("Finished running nlsLM #3\n")
            ##SIRX with automatic minimization to find kIX
            modelSIRXerr = class(modelSIRX) == 'try-error'
            if (!modelSIRXerr)
            {
                ##rm(beta0)
                coefsX2 = modelSIRX$m$getPars()
                start=c(X0=coefsX2[1], dX0 = coefsX2[2], beta=coefsX2[3], kIX=kIX, betaChgFact=betaChgFact)
                lower=c(X0=0, dX0=0,beta = 0, kIX=0.000001, betaChgFact=0.00)
                upper=c(X0=+N, dX0=0.5*N, beta=Inf, kIX=0.999, betaChgFact=1.0)
                
                names(start) = names(lower) = names(upper) = c('X0', 'dX0', 'beta', 'kIX', 'betaChgFact')
                 if (tChg<=0) 
                {
                    tChg = 0
                    betaChgFact = 1
                }
                if (!betaAuto)
                    {
                        start = start[c(1:4)]
                        lower = lower[c(1:4)]
                        upper = upper[c(1:4)]
                        ##betaChgFact = 1
                    }

                ##betaChgFact00 = betaChgFact
                

                ##browser()
                if (bAuto) printf("running modelSIRXauto\n")
                if (bTrace && bAuto) printf("Running nlsLM #4\n")
                if (bAuto) modelSIRXauto =
                               try(nlsLM(log(C*isC + dC*isDC)~ log(sirX3(
                                             retX=isC,
                                             retDX=isDC,
                                             tABtest=dTABcal,
                                             ABcurve=ABcurve,
                                             Tmax=dTmax,
                                             N=N,
                                             errorRetVal=(isC*C + isDC*dC),
                                             tChg = tChg,
                                             gamma=gamma,  ##gamma and above are known
                                             ##betaChgFact = betaChgFact0,
                                             ##X0=X0, dX0=dX0, beta=beta, beta0=beta0, kIX=kIX
                                             X0=X0, dX0=dX0, beta=beta, kIX=kIX, betaChgFact=betaChgFact
                                             ##,beta0=beta0, tChg=tChg  ##this is being solved for
                                         )),
                                         data=six2,
                                         start=start,
                                         lower=lower,
                                         upper=upper,
                                         weights = c(ABweight, rep(1, dim(six2)[1]-1)),
                                         algorithm='LM', trace=bTrace, control=list(warnOnly=T, maxiter=200, maxfev=500)
                                         ))
                if (bTrace && bAuto) printf("Finished nlsLM #4\n")

                ##browser()
                modelSIRXautoerr = class(modelSIRXauto) == 'try-error'
                ##THERE IS A BUG HERE - found it...tChg wasn't being set on the below sirX3 call
                if (!bAuto || modelSIRXautoerr)
                {
                    ##browser()
                    if (length(coefsX2) == 3)
                    {
                        if (bPrintOutOfContext)
                            printf("Drawing4param: tChg=%d X0=%f dX0=%f kIX=%f beta=%f betaChgFactIn=%f\n",
                                   tChg, coefsX2[1], coefsX2[2], kIX, coefsX2[3], betaChgFact)
                        xdat15rk = sirX3(N=N, X0=coefsX2[1], dX0=coefsX2[2], kIX=kIX,
                                         gamma=gamma, beta=coefsX2[3], Tmax=nProject+tMax-tMin, tChg=tChg,
                                         betaChgFact=betaChgFact)
                    }
                    else
                    {
                        if (bPrintOutOfContext)
                            printf("Drawing5param: tChg=%d X0=%f dX0=%f kIX=%f beta=%f betaChgFact=%f\n",
                                   tChg, coefsX2[1], coefsX2[2], kIX, coefsX2[3], coefsX2[4])
                        xdat15rk = sirX3(N=N, X0=coefsX2[1], dX0=coefsX2[2], kIX=kIX,
                                         gamma=gamma, beta=coefsX2[3], betaChgFact=coefsX2[4], Tmax=nProject+tMax-tMin, tChg=tChg)
                    }
                }
                else
                {
                    coefsX2a = modelSIRXauto$m$getPars()
                    ##browser()
                    if (length(coefsX2a) == 4)
                    {
                        kIXauto=coefsX2a[4]
                        if (bPrintOutOfContext)
                            printf("Drawing4autoKIXparam: tChg=%d X0=%f dX0=%f kIX=%f beta=%f betaChgFactIn=%f\n",
                                   tChg, coefsX2a[1], coefsX2a[2], kIXauto, coefsX2a[3], betaChgFact)
                         xdat15rk = sirX3(N=N, X0=coefsX2a[1], dX0=coefsX2a[2], kIX=kIXauto,
                                         gamma=gamma, beta=coefsX2a[3], Tmax=nProject+tMax-tMin, tChg=tChg, betaChgFact=betaChgFact)
                    }
                    else
                    {
                        kIXauto=coefsX2a[4]
                        if (bPrintOutOfContext)
                            printf("Drawing5autoKIXparam: tChg=%d, X0=%f dX0=%f kIX=%f beta=%f betaChgFact=%f\n",
                                   tChg, coefsX2a[1], coefsX2a[2], kIXauto, coefsX2a[3], coefsX2a[5])
                        ##browser()
                        xdat15rk = sirX3(N=N, X0=coefsX2a[1], dX0=coefsX2a[2], kIX=kIXauto,
                                         gamma=gamma, beta=coefsX2a[3], betaChgFact=coefsX2a[5], Tmax=nProject+tMax-tMin, tChg=tChg)
                    }   
          
                }
            }
                  
            nMissing0 = dim(xdat)[1] - dim(xdat15nlm)[1]
            nMissing = nMissing0
            if (nMissing0<0) nMissing = 0
             ##xdat=cbind(xdat, Spct=rep(NA, dim(xdat)[1]))
            xdat=cbind(xdat, Ix2=rep(NA, dim(xdat)[1]))
            xdat=cbind(xdat, Cx2=rep(NA, dim(xdat)[1]))
            xdat=cbind(xdat, dCx2=rep(NA, dim(xdat)[1]))
            xdat=cbind(xdat, Spct2=rep(NA, dim(xdat)[1]))
            ##dCsir = c(NA, xdat15nlm$X[-1] - xdat15nlm$X[-length(xdat15nlm$X)])
            dCsir = xdat15nlm$dX
            ##browser()
            ##if (nMissing>0)
            infected=xdat[,1:2]
            infected=cbind(infected, I=rep(NA, dim(infected)[1]))
            infected=cbind(infected, dI=rep(NA, dim(infected)[1]))
            for (i in 1:length(xdat$t))
            {
                it = xdat$t[i]
                ixdat15 = which(it == xdat15nlm$t + tMin)
                if (!modelSIRXerr) ixdat2  = which(it ==  xdat15rk$t + tMin)
                if (length(ixdat15 == 1) &&
                    FALSE) ##comment out for the Isir model
                {
                    xdat$Isir[i] = xdat15nlm$I[ixdat15]
                    xdat$Csir[i] = xdat15nlm$X[ixdat15]
                    xdat$dCsir[i] = dCsir[ixdat15]
                    ##xdat$Spct[i] = xdat15nlm$S[ixdat15]/N
                    ##xdat$Spct = round(xdat$Spct,4)
                    xdat$Isir = round(xdat$Isir,1)
                    xdat$Csir = round(xdat$Csir,1)
                    xdat$dCsir = round(xdat$dCsir,1)
##                    xdat$pct2 = round(xdat$Spct, 2)
                }
                if (!modelSIRXerr && length(ixdat2 == 1))
                {
                    ##browser()
                    xdat$Ix2[i] = xdat15rk$I[ixdat2]
                    xdat$Cx2[i] = xdat15rk$X[ixdat2]
                    xdat$dCx2[i] = xdat15rk$dX[ixdat2]
                    infected$I[i] = N - xdat15rk$X[ixdat2] - xdat15rk$S[ixdat2]
                    xdat$Spct2[i] = xdat15rk$S[ixdat2]/N
                    xdat$Spct2= round(xdat$Spct2,4)
                }
            }
            ##if (!modelSIRXerr)
            ##    xdat$dCx2 = c(NA, xdat$Cx2[-1] - xdat$Cx2[-length(xdat$Cx2)])

            ##xdat$Spct = 100*xdat$Spct
            xdat$Spct2 = 100*xdat$Spct2
            
         }
     }   
    
    
    rownames(xdat) = 1:nSim
    iData = xdat$t>=tMin & xdat$t<=tMax
    yrange = c(1, 2*max(na.rm=T, xdat$C, xdat$Csim, xdat$Cnlm))
    if (!is.na(N)) yrange[2] = min(N, yrange[2])

    yText = 10^(0:10)
    iyText = yText>=yrange[1] & yText<=yrange[2]
    yText = yText[iyText]
    par(xaxs = 'i', yaxs='i',oma=c(2,2,2,2.5))
    xticks = 5*(floor(min(na.rm=T,xdat$t)/5):ceiling(max(na.rm=T, xdat$t)/5))
    xdat$C[xdat$C<=0]=NA
    with(xdat, plot(C~t, log='y', type='o', pch='o', ylim=yrange, yaxt='n',
                    ylab = NA, xlab=NA,
                    xaxt='n', lab=c(8,30,20), main=country))
    text(x=min(xdat$t, na.rm=T), y=yText, labels=yText, pos=2, xpd=T) ##label y axis by hand
    xticksEven = xticks[xticks %% 2 == 0]
    xticksEven = xticksEven[xticksEven<=max(na.rm=T, xdat$t)]
    text(y=1, x=xticksEven,
         sprintf("%d\n%s", xticksEven, format(xdat$Date[1] + xticksEven - xdat$t[1], "%m-%d")),
         pos=1, xpd=T) ##label y axis by hand

    mtext("O's Cases, X's Daily Change (Black=data Blue=Calibration)",
          side = 2, line = -1, outer = T,
           adj = NA, padj = NA, cex = 1.2, col = 'black')

    
    mtext("t (Days past 12/31/19) and Date",
          side = 1, line = -2, outer = T,
           adj = NA, padj = NA, cex = 1.2, col = 'black')

    mtext("% Susceptible Population (Red)", side = 4, line = -1.75, outer = T, at = 0.33,
           adj = NA, padj = NA, cex = 1.2, col = 'darkred')

    
    
 ##    labels = ytick, srt = 45, pos = 2, xpd = TRUE)

    ##v=5*(floor(min(na.rm=T,xdat$t)/5):ceiling(min(na.rm=T, xdat$t)/5))
    if (any(xticks == tMin)) xticks = xticks[-which(xticks == tMin)]
    if (any(xticks == tMax)) xticks = xticks[-which(xticks == tMax)]

    lightticks = rgb(0.5,0.5,0.5)
    darkticks = rgb(0.3,0.3,0.3)
    
    abline(col=lightticks, v= xticks )
    abline(col=lightticks, h=c(10^(0:10), 2*10^(0:9), 4*10^(0:9), 6*10^(0:9), 8*10^(0:9) ) )
    abline(col=darkticks, h=c(10^(0:10)) )
    ##draw region of interest and tChg
    abline(col=rgb(0,0,1), v = c(tMin, tMax))
    if (tChg != 0)  abline(col=rgb(0,0.5,0), v = tMin + tChg)
    if (!is.na(tABcal)) abline(col=rgb(0.25,0.25, 0.25), v = tABcal) 
           

    nticks =  round(log10(max(yrange)/min(yrange))*2)
    ##grid(col=1, ny = nticks, equilogs=F)
         
    if (bPlotExtra)
        { ##For the paper can do without the non-linear model -- make chart less busy
            with(xdat, points(Cnlm~t, type='o', pch='+', col='darkgreen'))
            with(xdat, points(dCnlm~t, type='o', pch='+', col='darkgreen'))
        }

    if (bPlotExtra) with(xdat[iData,], points(Csim~t,  type='l', pch=0, col='darkred'))
    if (bPlotExtra) with(xdat[xdat$t<=tMin,], points(Csim~t, type='l', pch=0, col='darkblue'))
    if (bPlotExtra) with(xdat[xdat$t>=tMax,], points(Csim~t, type='l', pch=0, col='darkblue'))
    
    

    with(xdat, points(dC~t, type='o', pch='x'))
    if (bPlotExtra) with(xdat[iData,], points(dCsim~t, type='l', pch='+', col='darkred'))
    if (bPlotExtra) with(xdat[xdat$t<=tMin,], points(dCsim~t, type='l', col='darkblue'))
    if (bPlotExtra) with(xdat[xdat$t>=tMax,], points(dCsim~t, type='l', col='darkblue'))

    if(!is.null(bYtotInterval))
    {
                interval = match(bYtotInterval, six$t)
                with(xdat[interval,], points(dCold~t, type='o', pch='x', col=rgb(0.4,0.4,0.2)))
                with(xdat[interval,], points(Cold~t, type='l', pch='o', col=rgb(0.4,0.4,0.2)))
    }
    
    ##browser()

##    if (havFwd && any(colnames(xdat) == 'Csir'))
    if (any(colnames(xdat) == 'Csir'))
    {
        
        ##with(xdat, lines(Cfwd~t, col='purple'))
        ##with(xdat, lines(dCfwd~t, col='purple'))

        ##browser()
        if (bPlotExtra) with(xdat, points(Csir~t, col='red', type='l',pch='S'))
        if (bPlotExtra) with(xdat, points(dCsir~t, col='red', type='l', pch='S'))
    }

    if (!modelSIRXerr)
        {
            with(xdat, points(Cx2~t, col=rgb(0,0,1), pch='x', cex=0.7, type='o'))
            with(xdat, points(dCx2~t, col=rgb(0,0,1), pch='x', cex=0.7, type='o'))
        }

    ##to stop excessive digit printing and useless scientific notation
    xdat$Csim  = round(xdat$Csim, 1)
    xdat$dCsim = round(xdat$dCsim, 1)
    xdat$Cnlm  = round(xdat$Cnlm, 1)
    xdat$dCnlm = round(xdat$dCnlm, 1)

    xdat$Csim[xdat$Csim < -99999] = NA
    xdat$dCsim[is.na(xdat$Csim)] = NA
    xdat$Cnlm[xdat$Cnlm < -99999] = NA
    xdat$dCnlm[is.na(xdat$Cnlm)] = NA

    ##browser()
    if (!is.null(xdat$Spct2))
        with(xdat, points(Spct2~t, type='l', col='darkred', lwd=2)) ##if available plot percent remaining susceptibles
 
    ##printf("nticks=%d\n", nticks)

    xdat$C = round(xdat$C)
    xdat$dC = round(xdat$dC)
    ##browser()
    datagaps = which(xdat$t[-1] - xdat$t[-nSim] > 1)
    datagaps = datagaps[xdat$t[datagaps]>=tMin]

    xdat=xdat[xdat$t<=(tMax+nProject),]

    if (!bPrintOutOfContext) ##Simplify the table for the reader
        {
            xdatSimp = xdat[,c('Date','t','C','dC','Ix2','Cx2','dCx2','Spct2')]
            xdatSimp$Ix2 = round(xdatSimp$Ix2, 1)
            xdatSimp$Cx2 = round(xdatSimp$Cx2, 1)
            xdatSimp$dCx2 = round(xdatSimp$dCx2, 1)
            colnames(xdatSimp) = c('Date','t','C','dC','I(Model)','C(Model)','dC(Model)','S%(Model)')
            print(xdatSimp)
        }
    else
        print(digits=5, xdat)
    


    QprobAuto = kIXauto/(kIXauto + gamma) ##a doublecheck

    printf("-------------\n")
    printf("Country/State: %s\tData Source: %s\n", country, szDataSource[dataSource])
    if (bPrintOutOfContext)
        {
            printf("N=%.4g Aref=%g  Cdc=%g  rDC=%g To=%f Tmin=%d gamma=%f betaChgFact=%f\n", N, Aref, Cdc, rDC, Tref, tMin, gamma, betaChgFact)
    ##if (havFwd)
            printf("Reff=%f RinfPct=%f Xinf=%7.0f Io=%f kIX=%f Qprob=%f kIXauto=%g QprobAuto=%f\n",
                   Reff, RinfPct, Xinf, Io, kIX, Qprob, kIXauto, QprobAuto)
        }
     
    printf("---------------Log-Linear Bootstrap Model DC----------------\n")
    mySummary(modelDC)
    printf("---------------Log-Linear Bootstrap Model DC3---------------\n")
    mySummary(modelDC3)
    printf("---------------Closed Form Solution Exponential Model-------\n")
    mySummary(modelNls)
    printf("---------------SIRX Model Boostrap W/Preset Qprob-----------\n")
    mySummary(modelSIRX)
    if (bAuto)
        {
            printf("---------------SIRX Model Calibration Run-----------\n")  
            mySummary(modelSIRXauto)
        }
    if (!is.null(bYtotInterval)) mySummary(modYtot)

    if (any(datagaps))
        printf("Note! -- data gap causes non-linear dC after t=%d\n",
               xdat$t[datagaps] )


    if (!bPrintOutOfContext) ##Show the _final_ results -- make easy for readers
    {
        printf("------------------------\n")
        printf("Calibration Run Summary:\n")
        printf("Location=%s N=%.4g tMin=%d tMax=%d tBetaChange=%d tSeroTest=%d CasesDataSource=%s\n",
               country, N, tMin, tMax, tChg + tMin, tABcal, szDataSource[dataSource])

        coefsFinal=modelSIRXauto$m$getPars()
        printf("Gamma=%f Beta0=%f BetaChgFactor=%f Beta1=%f kIX=%f Qprob=%f\n",
               gamma, coefsFinal['beta'], coefsFinal['betaChgFact'],
               coefsFinal['beta']*coefsFinal['betaChgFact'],
               coefsFinal['kIX'],  coefsFinal['kIX']/( coefsFinal['kIX'] + gamma)
               )
    }
    if (!is.na(tABcal))
    {
        infected$I[is.na(infected$I)] = 0 
        infected$dI = c(NA, infected$I[-1] - infected$I[-dim(infected)[1]])
        i = which(infected$t<=tABcal)
        if (bTrace) print(infected[i,])
        NtestPop = N - xdat$Cx2[tABcal == xdat$t] ##This is correction for N as per Eq 12 of paper
        predictedTest = testPctFromNumInfected(infected$dI[i], tABcal - infected$t[i], nPop = NtestPop , curve=ABcurve)
        perfectTest = testPctFromNumInfected(infected$dI[i], tABcal - infected$t[i], nPop = NtestPop, curve=perfectAB)
        printf("Perfect/Predicted Test Results = %f / %f vs Calibration Target=%f\n", perfectTest, predictedTest, ABtarg)                
    }

    if (bPrintOutOfContext)
        printf("WARNING: PRINTING OUT OF CONTEXT EXTRA DEBUGGING INFORMATION ON -- READ THE CODE TO UNDERSTAND!!!")
    ##return(NULL)
    if (bRetxdat) return(xdat)
}




#############################################################################################
## SIR-X2 model after Maier and Brockmann https://doi.org/10.1101/2020.02.18.20024414 with certain changes
## Descretized in integer time steps using simple linear approximation. (a Runge-Kutta method would converge faster)
## The names variables are defined as in Maier et al.
## Note: in Maier et al the alpha is https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology beta.
## Maier's beta is wikipedia's gamma.
## (c) Michael Halem 2020  BSD License (i.e. Do what you want with it, but please attribute)
## If you use the code I would appreciate a small citation in your work in addition to Ben Maier's (who invented the concept).
## The below initial values are from Maier et al's ex-Hubei. I have a few suggestions on modifications to fit better the
## underlying physics (i.e. reality) of the epidemic.  You may email me at my public email address: consult@rocketpwr.com
## t = Days since start of simulation
## S = Succeptible
## I = Infected
## X = Cumulative Reported Cases
## Additional outputs are:
## RfromS = Removed directly from susceptible compartment without infection
## RfromI = Removed from Infected Compartment without passing thru Confirmed Cases
## This is essentially a sorting of both infected and susceptibles into impervious X and RfromS compartments.
## An improvement would be to make the RfromS not impermeable, but allow a different, lower number of infections from
## that compartment so that d[RfromS] = ... -delta * SrFromS[i0]/N*I[i0] and
## d[I] = + delta *  SrFromS[i0]/N*I[i0]
## this delta in the same units and analogous to beta,
## renaming for the SPIXR (spiker) as (SPICR is offensive) where new compartment P is susceptibles taking precautions
## P has a new beta equavalent, delta, which is the amount of people who get infected from one infected in their compartment
## Note: model still assumes perfect quarentine of confirmed cases

sirX2<-function(N=1e6, X0=70,
               I0.X0=6.98, ##Note: A=1 => Io/Xo = 0 means 100% are quarantined => growth = 0
               kSP= 0,  ## S to Protected
               kIX= 0.125, ##jump from X to the quarantined pen
               gamma = 0.125,  ##my conventional beta and gamma, i.e. gamma is from I to R compartment
               beta = 0.375,  ##from S to I compartment
               delta = 1/2 * gamma, ##from P (precautionarys) to infected 
               Tmax=45,
               A = NA, ##Ascertainment ratio, i.e. number of Actual infected to reported infected = (I + X)/X = 1 + X/I
               dt=0.001,
               bRound=T) ##round results
{

    nObs = round(Tmax/dt)
    nStep = round(1/dt)
    
    Ro.free = beta/gamma
   
    
    
    Pben = kSP/kIX
    Q = kIX/(gamma + kIX)

    t = 0
    P = 0 ##protecteds / paranoids
    ##RfromI = (beta/(k0+k) - 1) * X0 ##removed via infection
    R = 0 ## This is an unobserved integration constant and doesn't matter in long run
    X = X0    ##confirmed cases
    if (!is.na(A)) I0.X0 = A - 1
    I = X0 * I0.X0 ##infections
    S = N - P - R - X - I
    
    ##descrete time steps with no correction for integer step sizes
    
    for (i in 2:(nObs+1))
    {
        i0 = i - 1
        S[i] = S[i0] + (- beta * S[i0]/N * I[i0] - kSP * S[i0]) * dt
        I[i] = I[i0] + (beta *S[i0]*I[i0]/N - gamma*I[i0] - (kIX)*I[i0]) * dt + delta * P[i0]/N * I[i0] * dt
        X[i] = X[i0] + ((kIX) * I[i0]) * dt
        P[i] = P[i0] + (kSP * S[i0]) * dt - delta * P[i0]/N * I[i0] * dt
        R[i] = R[i0] + (gamma * I[i0]) * dt
        t[i] = (i -1) * dt
    }

    printf("SIX-X N      = %10.0f X0    = %7.0f Io/Xo = %7.4f\n", N, X0, I0.X0)
    printf("      Pben   = %10.4f Q     = %7.4f A     = %7.4f\n", Pben, Q, A) 
    printf("      Ro.free= %10.4f beta  = %7.4f gamma = %7.4f\n", Ro.free, beta, gamma)
    printf("      R.prot = %10.4f delta = %7.4f\n", delta/gamma, delta)    
    printf("      kSP    = %10.4f kIX   = %7.4f\n", kSP, kIX)
    printf("------------\n")

    
    
    dat = data.frame(t, S, P, I, X, R)
    if (bRound) dat = round(dat,0)
    dat = dat[(0:Tmax)*nStep + 1,]
    return(dat)
}



#####################################################################
## testPctFromNumInfected
## Given the true number of seropositive time series and the population
## return the percentage that would be tested positive by an antibody test
## with a known input curve

testPctFromNumInfected<-function(nSeroPos, t, nPop = NA, curve=perfectAB)
{
    nTestPos = 0
    nFalsePos = 0

    if (length(t) <= 0) return(NA)
    if (length(t) != length(nSeroPos)) return(NA)

    nSeroPos[is.na(nSeroPos)] = 0
    for (j in 1:length(t))
        {
            i = max(which(t[j] >= curve$t0)) ##Note: can fit a curve, but beyond scope of current analysis, i.e. p = 1 - exp(r*t) curve
            p = curve$sens[i]
            q = curve$spec[i]
            nTestPos = nTestPos + nSeroPos[j] * (p + q - 1)
            nFalsePos = nTestPos * (1 - q)  ##just making a weighted average of some kind
        }
    
    ##browser()
    
    pctPos = nTestPos/nPop + nFalsePos/nTestPos
            
    return(pctPos)


}

######################################################################
## Make the antibody surve plots for the paper

makeAntiBodyCurve<-function(abcurve=exampleAB)
{
    N = 10000 ##Note R Language doesn't like the naked N as is shared with the next command
    t0 = 40
    t = 0:t0
    dt = t0 - t
    r = log(2)/3
    Iexp = pmin(N,floor(exp(t*r)))
    dIexp = c(Iexp[1], Iexp[-1] - Iexp[-(t0+1)])
    dIdecay = dIexp[(t0+1):1]
    Idecay = rep(NA, t0+1)
    Idecay[1] = dIdecay[1]
    for (i in 2:(t0+1)) Idecay[i] = Idecay[i-1]+dIdecay[i]
    Iflat = dIflat = rep(NA, t0+1)
    dIflat[1:(t0+1)] = ceiling(N/(t0+1))
    ##browser()
    dIflat[t0+1] = N - sum(dIflat[1:t0])
    Iflat[1] = dIflat[1]
    for (i in 2:(t0+1)) Iflat[i] = Iflat[i-1]+dIflat[i]

    p = approx(x=abcurve$t0, y=abcurve$sens, xout=dt, method='constant', rule=2)$y
    q = approx(x=abcurve$t0, y=abcurve$spec, xout=dt, method='constant', rule=2)$y
    q = mean(q)
    
    
    TleftExp = (p + q - 1)*dIexp
    TleftFlat = (p + q - 1)*dIflat
    TleftDecay = (p + q - 1)*dIdecay
    

    results = cbind(t,dt,Iexp, dIexp, TleftExp, Iflat, dIflat, TleftFlat, Idecay, dIdecay, TleftDecay)
    rownames(results) = 1:dim(results)[1]

    thetaExp = sum(TleftExp)/N + 1 - q
    thetaFlat = sum(TleftFlat)/N + 1 - q
    thetaDecay = sum(TleftDecay)/N + 1 - q

    results = as.data.frame(results, stringsAsFactors=F)
    print(results)

    printf("thetaExp=%.3f%%, thetaFlat=%.3f%%, thetaDecay=%.3f%%\n", 100*thetaExp, 100*thetaFlat, 100*thetaDecay)

    ##with(results, plot(y=dIexp, x=t, log='y',ylim=c(yMin,150), xlim=c(0,40), col='black', type='S'))

    
    plotIt(y=results$Iexp, dy=results$dIexp, t=results$t, N=N, theta=thetaExp, filename="expPlot.pdf",
           name="Figure 1 - Exponential Increase of New Infected Time Series")
    dev.new()
    
    plotIt(y=results$Iflat, dy=results$dIflat, t=results$t, N=N, theta=thetaFlat, filename="flatPlot.pdf",
           name="Figure 2 - Flat New Infected Time Series")
    
    dev.new()
    
    plotIt(y=results$Idecay, dy=results$dIdecay, t=results$t, N=N, theta=thetaDecay, filename="decayPlot.pdf",
           name="Figure 3 - Exponential Decrease of New Infected Time Series")

    
}
    
####################################################################################################
## Plot barchart and change chart for paper


plotIt<-function(y, dy, t, N, theta, name="???", filename="theFile", dir="~/Downloads")
    {
        yMin = 0.0080
        yMax = 150
        yTicks = c(0.01, 0.1, 1, 10, 100)

        ##convert to percent
        y = 100 * y/N
        dy = 100 * dy/N
        pdf(sprintf("%s/%s", dir, filename))
        
        plot(1, type='n', xlab='t (Time in Days)', main=name,
             ylab="X's = New True Positives, Bars = Cumulative True Positives ",
             log='y',ylim=c(yMin,150), xlim=c(0,40), yaxt='n', xaxs='i', yaxs='i')
        axis(2, at=yTicks, lab=sprintf("%.2f%%",yTicks))
        for (i in 2:length(y))
            polygon(y=c(yMin,yMin,y[i],y[i],yMin), x=c(t[i-1], t[i], t[i], t[i-1], t[i-1]), col='white' )
        points(y=dy, x=t, col='black', type='o', pch='X')
        abline(h=100*theta, col='darkred')
        abline(h=100, col='darkblue')
        graphics.off()                                                                  
     }

