library(glmnet)
library(Matrix)
library(MASS)
#rm(list=ls())
#Rcpp::sourceCpp('C:/Users/Knight/Dropbox/yz_lq/new_0113.cpp')
Rcpp::sourceCpp("E:/ZeroinflatedCode/R/newest_Newton_finalvison.cpp")
set.seed(2022)


TSIA<-function(p,n,n1,m,threshold,printyes,usealpha,usetest,repeatnumber,setting,betareal,H){
  twostage <- function(X, Y, Xtest, Ytest){
    indexnonzero = which(Y != 0)
    Xnew = as.matrix(t(X[,indexnonzero]))
    logY = as.matrix(log(sqrt(Y[indexnonzero]/(1-Y[indexnonzero]))),1,length(indexnonzero))
    BlogY = as.matrix(Y[indexnonzero],1,length(indexnonzero))
    solution = glmnet(Xnew, t(logY) , family="gaussian")

    SYhat=predict(solution,t(X),s=0.00000,type ="response")

    Yhat<-(exp(SYhat))^2/(1+(exp(SYhat))^2)

    Ytrain = Y
    Ytrain[Y != 0] = 1
    glmmod=glmnet(t(X),t(Ytrain),family="binomial")
    classify<-predict(glmmod,t(X),s=0.00000,type='response')
    jkl<-glmmod[["beta"]]
    jkll<-as.matrix(jkl)
    alpha<-jkll[,ncol(jkll)]
    class<-ifelse(classify<=0.5,0,1)
    indexzero = which(class == 0)
    Yhat[indexzero] = 0

    objtrain = norm(Y-Yhat, "2")^2 / length(Y)



    SYhattest=predict(solution,t(Xtest),s=0.00000,type ="response")

    Yhattest<-(exp(SYhattest))^2/(1+(exp(SYhattest))^2)
    classifytest<-predict(glmmod,t(Xtest),s=0.00000,type='response')
    classtest<-ifelse(classifytest<=0.5,0,1)
    indexzerotest = which(classtest == 0)
    Yhattest[indexzerotest] = 0

    objtest = norm(Ytest-Yhattest, "2")^2 / length(Ytest)


    return(list(objtrain =objtrain, objtest = objtest,beta = solution$beta,alpha=alpha))
  }
#p = 6;
#n = 400;
#n1 = 100;
#m = 50;
#threshold = -1;
#printyes = TRUE
#usealpha = FALSE
#usetest = TRUE
#repeatnumber = 10
#setting = 1
#betareal = 0.2*matrix(c(-2, 2,-2.5, 1, -1, 2),p,1 )
#betareal = 0.3*matrix(c(-3.791075, 4.009764, -4.767333, 1.813433, -2.448710, 4.545507), p,1)
#betareal = 0.3*matrix(c(-3, 5, -2, 3, -1, 4), p,1)
#betareal = 0.1*matrix(c(-3, 5), p,1)
#betareal = matrix(rnorm(p), p,1)
#if (betareal[1] > 0){
#  betareal[1] = - betareal[1]
#}

cat("\n realbeta = ", betareal)
cat("\n p = ", p, ", n = ", n, ", m = ", m)
cat("\n threshold = ", threshold)
cat("\n repeatnum = ", repeatnumber)


sample_cov = matrix(rnorm((p+2)*p),(p+2),p)
test_cov = cov(sample_cov)



#H = c(0.1)
#H = 0.1
#H = seq(from = 0.01, to = 0.09, by = 0.01)
finalobj = matrix(0,length(H),1)
finalobj2 = matrix(0,length(H),1)
finalobj3 = matrix(0,length(H),1)

# model 1
averagetest = matrix(0,p,length(H))
averageteststep1 = matrix(0,p,length(H))
averageobj = matrix(0,length(H),1)
averageobjtest = matrix(0,length(H),1)

# model 2
averagetest_model2 = matrix(0,p,length(H))
averageobj_model2 = matrix(0,length(H),1)
averageobjtest_model2 = matrix(0,length(H),1)


averageobjteststep1 = matrix(0,length(H),1)

# two stage
averagetest2 = matrix(0,p,1)
averageobj2 = 0
averageobj2test = 0

averageobj2testB = 0
#for (k in 1:length(H)) {
#  h = H[k]
#  cat("\n\n *********************************\n  h = ", h, "\n")

for(i in 1:repeatnumber){
  Xall = matrix(rnorm(p*(n+n1), mean = 1), p, n+n1);
  #Xall = matrix(rbeta(p*(n+n1), 3,2 ), p, n+n1);
  #Xall[1,] = 1
  XTbetaall = t(Xall) %*% betareal;

  if (setting == 1){
    prob = (exp(XTbetaall))^2 / (1 + (exp(XTbetaall))^2)
    Yall = matrix(rbinom(n+n1,m,prob),n+n1,1) / (m+1)
  }else if (setting == 2){
    prob = 1 - exp(- XTbetaall)
    index_great_prob = which(prob > 0)
    great_number = length(index_great_prob)
    Y_great = matrix(rbinom(great_number,m,prob[index_great_prob]), great_number, 1) / (m + 1)
    Yall = matrix(0, n+n1,1)
    Yall[index_great_prob] = Y_great
  }
  Yall[which(XTbetaall <= threshold)] = 0

  X = Xall[,1:n]
  Y = Yall[1:n,1]

  Xtest = Xall[,(n+1):(n+n1)]
  Ytest = Yall[(n+1):(n+n1),]



  solutiontwostage = twostage(X,Y,Xtest,Ytest)
  betatwostage = as.matrix(solutiontwostage$beta)
  betatwostage = betatwostage[,ncol(betatwostage)]
  averagetest2 = averagetest2 + betatwostage
  averageobj2 = averageobj2 + solutiontwostage$objtrain
  averageobj2test = averageobj2test + solutiontwostage$objtest


  averageobj2testB=averageobj2testB+solutiontwostage$Bobjtest
  alphatwo = as.matrix(solutiontwostage$alpha)



  stoptol = 1e-7
  maxiter = 50
  beta = matrix(rnorm(p),p,1)
  beta1 = matrix(1,p,1)
  # alpha should be a p \times 1 matrix

  alpha = matrix(1,p,1)
  #alpha = alphatwo

  for(j in 1:length(H)){

    h = H[j]
    test = mytestnewest(X,Y,beta=beta1,h,stoptol,maxiter,printyes,alpha,usealpha,Xtest,Ytest,usetest)

    averagetest[,j] = averagetest[,j] + test$beta
    averageobj[j,1] = averageobj[j,1] + test$objtrain
    averageobjtest[j,1] = averageobjtest[j,1] + test$objtest
    averageobjteststep1[j,1] = averageobjteststep1[j,1] + test$objteststep1
    averageteststep1[,j] = averageteststep1[,j] + test$betastep1
  }



}

# model 1
averagetest = averagetest / repeatnumber
averageteststep1 = averageteststep1 / repeatnumber
averageobj = averageobj / repeatnumber
averageobjtest = averageobjtest / repeatnumber
averageobjteststep1 = averageobjteststep1 / repeatnumber

# model 2
averagetest_model2 = averagetest_model2 / repeatnumber
averageobj_model2 = averageobj_model2 / repeatnumber
averageobjtest_model2 = averageobjtest_model2 / repeatnumber


# two stage
averagetest2 = averagetest2 / repeatnumber
averageobj2 = averageobj2 / repeatnumber
averageobj2test = averageobj2test / repeatnumber


for (i in 1:length(H)) {
  cat("\n\n *********************** \n  h = ", H[i])
  cat("\n betareal = ", betareal, "\n")
  cat("\n   our method: \n     beta = ", averagetest[,i], "\n     objtrain =  ", averageobj[i,1], "\n     objtest =  ", averageobjtest[i,1])
  cat("\n   step 1: \n     beta = ", averageteststep1[,i], "\n     objtest =  ", averageobjteststep1[i,1])
  cat("\n   two stage: \n     beta = ", averagetest2, "\n     objtrain =  ", averageobj2, "\n     objtest =  ", averageobj2test )

}
}
#}

