###############################################################################################
#ChIPmix.R 07-09-25 Aim : Application of linear regression mixture
#model to ChIP-chip data analysis
#-----------------------------------------------------------------------------------------------
#Main function : - initialisation of the mixreg function with a PCA or randomly
#                - regression with one component
#                - regression with two  components (same variances)
#                - best model selected with BIC criterion
#                - choice of the enriched population using the intersection of the two
#                  regressions
#                - threshold and posterior probability definitions
#                  for each probe (False discovery control)
#-----------------------------------------------------------------------------------------------
#Arguments :     - fileIN = name of the input data file with 3 columns
#                  named ID, INPUT et IP (ID=identifier of the probe, INPUT =
#                  log2(INPUT), IP = log2(IP)). 
#                - eps.mixreg = eps parameter of mixreg function specifying
#                  the convergence criterion for the EM algorithm.   
#                - alpha = level to control false discovery rate
#                - init.by.PCA : mixreg initialisation. If TRUE
#                  mixreg initialisation is done with a PCA and randomly otherwise
#                  (default value = TRUE)
#                - fileOUTscreen = name of the text file for results (parameters estimation...)  #                - fileOUTgraph = name of the plot file
#                - fileOUT = name of the output data file with 7 columns
#                  (ID=probe identifier, INPUT = log2(INPUT), IP = log2(IP),
#                  tau = posterior   probability, threshold=s, status = probe status,
#                  1 if the probe is classified as enriched, 0 otherwise
#------------------------------------------------------------------------------------------------#
#Output :        - 3 files: one for the parameters estimation of the mixreg function,
#                           one for   the graph with the enriched probes colored in red
#                           and the two regression lines, and the last one for the analysis
#                           summary of all
#------------------------------------------------------------------------------------------------#
#Dependancies : library mixreg (package mixreg), use of the function "mixreg"
#------------------------------------------------------------------------------------------------#
#Authors : Caroline Bérard: caroline.berard@agroparistech.fr
#          Marie-Laure Martin-Magniette: marie_laure.martin@agroparistech.fr
#----------------------------------------------------------------------------------------------
#
#Date : 25/09/2007
#-----------------------------------------------------------------------------------------------
#License : GPL
#-----------------------------------------------------------------------------------------------
#Warning : The name of the files "fileIN" and "fileOUTscreen" have to
#          be written completely (with extension, .txt for example),
#          but not the "fileOUTgraph" and "fileOUT" files.
################################################################################################
ChIPmix = function(fileIN="data.txt",eps.mixreg=1e-06,alpha=0.01,init.by.PCA=TRUE,
  fileOUTscreen="Regression_Results.txt",fileOUTgraph="Graph",fileOUT="Output")
  {
    library(mixreg)
    
    data = read.table(fileIN,h=T);
    ID = data$ID;
    INPUT = data$INPUT;
    IP = data$IP;
    n=length(INPUT);
    x=range(INPUT);
    
    if (init.by.PCA == T)
      {
        #Initialisation for mixreg function : PCA 
        M = cbind(INPUT,IP);
        elempropres = eigen(cov(M));
        #y = a + bx
        b = elempropres$vectors[1,1]/elempropres$vectors[2,1];
        a = mean(IP) - b*mean(INPUT);
        var = 0.5;

        eps1 = 5*b/100;
        b1 = b + eps1;
        b2 = b - eps1;

        TS1 =  list(list(beta=c(a,b),sigsq=var,lambda=1));
        TS2 =  list(list(beta=c(a,b1),sigsq=var,lambda=1/2),
          list(beta=c(a,b2),sigsq=var,lambda=1/2));
      } else {TS1 = NULL; TS2 = NULL;}
    
#linear regression with 1 component in the mixture fit1 =
    mixreg(INPUT,IP,theta.start=TS1,ncomp=1,eps=eps.mixreg,eq.var=T,verb=F);

#linear regression with 2 components in the mixture
    fit2 = mixreg(INPUT,IP,ncomp=2,theta.start=TS2,eps=eps.mixreg,eq.var=T,verb=F);

    sink(fileOUTscreen);

    bic1 = -2*fit1$log.like + log(n)*3; 
    print("bic1 : ");
    print(bic1);
    cat("\n");  
    
    bic2 = -2*fit2$log.like + log(n)*6; 
    print("bic2 : ");
    print(bic2);
    cat("\n");
    
#Selected model with only one population 
    if(bic2 > bic1)
      {
                                        #Graph 
        pdf(file=paste(fileOUTgraph,".pdf",sep=""));
        plot(INPUT,IP,pch='.',xlab="INPUT",ylab="IP");
        lines(x,y=fit1$parmat[1,1] + fit1$parmat[1,2]*x,col="red")
        dev.off()
        
        print("EM convergence :");
        print(fit1$converged);
        cat("\n");
        print("nb step :");
        print(fit1$nsteps);
        cat("\n");
        print("log.like :");
        print(fit1$log.like);
        cat("\n");
        parmat = matrix(t(fit1$parmat),byrow=T,ncol=4,nrow=1
          ,dimnames=list(c("Group0"),c("a","b","std.error","pi")))
        print(parmat);
        cat("\n");
        invisible(parmat)
      }
#Selected model with two populations     
    if (bic2 <= bic1)
      {
    #definition of the enriched population 
        fit2$parmat = fit2$parmat[order(fit2$parmat[,4],decreasing=T), ];
        
        a0 = fit2$parmat[1,1];
        a1 = fit2$parmat[2,1];
        b0= fit2$parmat[1,2];
        b1 = fit2$parmat[2,2];
    #intersection coordinates of the 2 lines
        abs = (a0-a1)/(b1-b0);
        probe.number.larger.abs = sum(INPUT>abs);
        probe.number.smaller.abs = n - probe.number.larger.abs;
        if (probe.number.larger.abs > probe.number.smaller.abs )
          {
            y1 = a1 + b1*(abs+1);
            y0 = a0 + b0*(abs+1);
            
          } else {
            y1 = a1 + b1*(abs-1);
            y0 = a0 + b0*(abs-1);
          }
     #the parameters of the enriched population are represented in
     #   the second row of the matrix "fit2$parmat"
        if (y0 > y1)
          {
            fit2$parmat = fit2$parmat[order(fit2$parmat[,4],decreasing=F), ];
            a0 = fit2$parmat[1,1];
            a1 = fit2$parmat[2,1];
            b0= fit2$parmat[1,2];
            b1 = fit2$parmat[2,2];
          }
        std.error0 = sqrt(fit2$parmat[1,3]); 
        std.error1 = sqrt(fit2$parmat[2,3]);
        pi = fit2$parmat[2,4];

    #definition of the posterior probability for a given probe to be enriched
        phi0 = dnorm(IP,a0+b0*INPUT,std.error0);
        phi1 = dnorm(IP,a1+b1*INPUT,std.error1);
        deno = (1-pi)*phi0 + pi*phi1;
        tau = pi*phi1 / deno;
                   
    #False discovery control : s is the threshold, which depends on INPUT and alpha.
        mu0 = (a0 + b0*INPUT);
        mu1 = (a1 + b1*INPUT);
        delta = ((mu1-mu0)/std.error0);
        lambda = delta*(qnorm((1-alpha),0,1)-delta/2) - log((1-pi)/pi);
        s = exp(lambda) / (1 + exp(lambda));
        if (probe.number.larger.abs > probe.number.smaller.abs)
          {
         #ind=1 if tau>s, 0 if not
            ind=0*numeric(length(tau));
            index<-which(tau-s>=0 & INPUT>abs)
            ind[index]<-1
          } else {
         #ind=1 if tau>s, 0 if not
            ind=0*numeric(length(tau));
            index<-which(tau-s>=0 & INPUT<abs)
            ind[index]<-1
          }


    #Color in red the enriched probes (ind=1) and the others in black (ind=0)
        couleur=c();
        couleur[ind == 1] = "red";
        couleur[ind == 0] = "black";
        
    #Output file
        S = data.frame(ID,INPUT,IP,tau,threshold=s,status=ind);
        write.table(S,file=paste(fileOUT,".txt",sep=""),sep="\t",row.names=FALSE);
      
    #Graph
        pdf(file=paste(fileOUTgraph,".pdf",sep=""));
        plot(INPUT,IP,pch='.',xlab="INPUT",ylab="IP",col=couleur);
        lines(x,y=b0*x + a0,col="red")
        lines(x,y=b1*x + a1,col="blue");
        dev.off()

        print("EM convergence :");
        print(fit2$converged);
        cat("\n");
        print("nb step :");
        print(fit2$nsteps);
        cat("\n");
        print("log.like :");
        print(fit2$log.like);
        cat("\n");
        parmat = matrix(t(fit2$parmat),byrow=T,ncol=4,nrow=2,dimnames=list(c("Group0","Group1"),c("a","b","std.error","pi")));
        print(parmat);
        cat("\n");
        invisible(parmat);
      }
    sink()
  }




