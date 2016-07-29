require ("MASS");
require ("e1071")
require ("randomForest")
require ("glmnet")
 
##
## @fn fitSingle (x, a, y, surrogate, augmodel=list(aug=F, predictor = NULL),
##                psmodel = list(ps=F,predictor= NULL),  cv=list(cv = T, folds = 5,
##                lambdas=2^seq(-5,5,1)))
## @brief Minimizes the IPWE (or model-assisted IPWE, under augmentation)
##        to contstruct a linear-thresholding decision rule.
## @param x n x p design matrix DOES NOT INCLUDE INTERCEPT
## @param a txt vector of length n txts (coded -1, 1)
## @param y outcome vector of length n coded so that higher is better
## @param surrogate options include 'logit' 'exp' 'sqhinge' 'hinge'
## @param augmodel: list indicating if an augmentation term should be used, if aug = False (default) then the augmentation term equals 0. Otherwise if True and aug='LR' then a ridged linear regression is used to estimate the augmentation term with the specified predictor; if True and aug='ML' then a random forest is used to estimate the augmentation term with the specified predictor.
## @param psmodel: list indicating if a propensity score model should be fitted, if ps = False (default) then the propensity score equals to 0.5 for all subjects. Otherwise if True and ps='LogitR' then a penalized logistic regression is used to estimate the propensity score with the specified predictor; if True and ps='ML' then a random forest is used to estimate the propensity score with the specified predictor.
## @param cv: list indicating if cross validation should be used. If cv = True (defualt), then cross validation will be used to select the optimal tuning parameter.
## @return list with the following key-->value pairs
##         bTilde--> stimated optimal beta, length p+1
##         wpos1/wneg1 --> weights

fitSingle = function (x, a, y, surrogate, augmodel=list(aug=F, predictor = NULL), 
         psmodel = list(ps=F,predictor= NULL), cv=list(cv = T, folds = 5, 
         lambdas=2^seq(-5,5,1))){
    aug = augmodel$aug
    ps  = psmodel$ps
    psx = psmodel$predictor
    regx = augmodel$predictor
         
    n = dim(x)[1];
    p = dim(x)[2];
    regp = dim(regx)[2];
    
    ## Fit propensity model
    phat = 0.5
    if (ps=='LogitR'){
		cvpsL2 <- cv.glmnet(psx, (a+1)/2, family='binomial',alpha=0,lambda=2^seq(-5,5,1))
		psprobL2 <- predict(cvpsL2,psx,s='lambda.min',type='response')
		phat <- psprobL2*(a==1)+(1-psprobL2)*(a==-1)
    }

    if (ps=='ML'){
        rfPS=randomForest (cbind (1, psx), as.factor(a))
        psprob = rfPS$votes
        Aind = which(colnames(psprob)=='1')
        psprobRF = psprob[,Aind]
        phat <- psprobRF*(a==1)+(1-psprobRF)*(a==-1)
    } 
    
    offset1 = rep (0, n);
    offset2 = rep (0, n);    
    if (aug == 'LR'){
        lamSeq = seq (0.05, 100, by=.25);
        d = cbind (1, regx, a, a*regx);
        lmrEst = lm.ridge (y ~ d-1, lambda=lamSeq);
        lmrCoeff = coefficients (lmrEst)[which.min(lmrEst$GCV),];
        alph = lmrCoeff[1:(regp+1)];
        beta = lmrCoeff[(regp+2):(2*regp+2)];
        offset1 = cbind (1, regx)%*%alph +  cbind(1,1*regx)%*%beta;
        offset2 = cbind (1, regx)%*%alph +  cbind(-1,-1*regx)%*%beta;
    }
    
    #    if (aug == 'LogitR'){
    #    d = cbind (1, regx, a, a*regx);
    #    cvaug = cv.glmnet(d,y,family='binomial',alpha=0,lambda=2^seq(-5,5,1))
    #    offset1 = predict(cvaug,cbind (1, regx, 1, regx),s='lambda.min',type='response')
    #    offset2 = predict(cvaug,cbind (1, regx, -1, -regx),s='lambda.min',type='response')
    #}
    
    if (aug == 'ML'){
        rfMdl = randomForest (cbind (1, regx, a, a*regx), y);
        offset1 = predict (rfMdl, cbind (1, regx, 1, 1*regx));
        offset2 = predict (rfMdl, cbind (1, regx, -1, -1*regx));
    }
    
    ## Recode txts and outcome so that weights are always positive
    y1 = (a==1)*(y -offset1)/phat + offset1
    y2 = (a==-1)*(y -offset2)/phat + offset2
     
    ### different surrogate functions 

	## default lambda
	lambda=0.5;

	###
	if (cv$cv>0){
		folds = cv$folds
		lambdas = cv$lambdas
		value_lambda = rep(0,length(lambdas))
		samplen = sample(n);
		for (k in 1:length(lambdas)){ 
			lambda = lambdas[k];
			value=0;
			for (i in 1:folds){
				tst_idx = samplen[seq(i,n,by=folds)]
				trn_idx = setdiff(1:n,tst_idx)
				objFn = function (b){
					if (surrogate == 'logit'){
						return (sum(y1[trn_idx]*(y1[trn_idx]>0)*
						log(1+exp(-b[1]-(x[trn_idx,]%*%
						b[-1]))))+sum(-y1[trn_idx]*(y1[trn_idx]<0)*
						log(1+exp(b[1]+(x[trn_idx,]%*%
						b[-1]))))+sum(y2[trn_idx]*(y2[trn_idx]>0)*
						log(1+exp(b[1]+(x[trn_idx,]%*%
						b[-1]))))+sum(-y2[trn_idx]*(y2[trn_idx]<0)*
						log(1+exp(-b[1]-(x[trn_idx,]%*%
						b[-1]))))+lambda*sum(b[-1]^2));
						}
					if (surrogate == 'exp'){
						return (sum(y1[trn_idx]*(y1[trn_idx]>0)*
						exp(-b[1]-(x[trn_idx,]%*%
						b[-1])))+sum(-y1[trn_idx]*(y1[trn_idx]<0)*
						exp(b[1]+(x[trn_idx,]%*%
						b[-1])))+sum(y2[trn_idx]*(y2[trn_idx]>0)*
						exp(b[1]+(x[trn_idx,]%*%
						b[-1])))+sum(-y2[trn_idx]*(y2[trn_idx]<0)*
						exp(-b[1]-(x[trn_idx,]%*%
						b[-1])))+lambda*sum(b[-1]^2));
						}	
					if (surrogate == 'sqhinge'){
						return (sum(y1[trn_idx]*(y1[trn_idx]>0)*
						(pmax(1-b[1]-(x[trn_idx,]%*%
						b[-1]),0))^2)+sum(-y1[trn_idx]*(y1[trn_idx]<0)*
						(pmax(1+b[1]+(x[trn_idx,]%*%
						b[-1]),0))^2)+sum(y2[trn_idx]*(y2[trn_idx]>0)*
						(pmax(1+b[1]+(x[trn_idx,]%*%
						b[-1]),0))^2)+sum(-y2[trn_idx]*(y2[trn_idx]<0)*
						(pmax(1-b[1]-(x[trn_idx,]%*%
						b[-1]),0))^2)+lambda*sum(b[-1]^2));
						}					
					if (surrogate == 'hinge'){
						return (sum(y1[trn_idx]*(y1[trn_idx]>0)*
						(pmax(1-b[1]-(x[trn_idx,]%*%
						b[-1]),0)))+sum(-y1[trn_idx]*(y1[trn_idx]<0)*
						(pmax(1+b[1]+(x[trn_idx,]%*%
						b[-1]),0)))+sum(y2[trn_idx]*(y2[trn_idx]>0)*
						(pmax(1+b[1]+(x[trn_idx,]%*%
						b[-1]),0)))+sum(-y2[trn_idx]*(y2[trn_idx]<0)*
						(pmax(1-b[1]-(x[trn_idx,]%*%
						b[-1]),0)))+lambda*sum(b[-1]^2));
						}						
				}
				#bTilde = optim (rnorm(p+1), objFn, method="BFGS")$par;
				bTilde = optim (rep(1,p+1), objFn, method="BFGS")$par;
				pred_tst = sign(cbind(1,x[tst_idx,])%*%bTilde)
				value = value+sum(y[tst_idx]*(a[tst_idx]==pred_tst))/sum(a[tst_idx]==pred_tst)
				value_lambda[k]=value
			}
		}
		lambda=lambdas[which.max(value_lambda)];
	}

	objFn = function (b){
		  if (surrogate == 'logit'){  
			return (sum(y1*(y1>0)*log(1+exp(-b[1]-(x%*%b[-1]))))+
			sum(-y1*(y1<0)*log(1+exp(b[1]+(x%*%b[-1]))))+
			sum(y2*(y2>0)*log(1+exp(b[1]+(x%*%b[-1]))))+
			sum(-y2*(y2<0)*log(1+exp(-b[1]-(x%*%b[-1]))))+lambda*sum(b[-1]^2))
			}
		  if (surrogate == 'exp'){  
			return (sum(y1*(y1>0)*exp(-b[1]-(x%*%b[-1])))+
			sum(-y1*(y1<0)*exp(b[1]+(x%*%b[-1])))+
			sum(y2*(y2>0)*exp(b[1]+(x%*%b[-1])))+
			sum(-y2*(y2<0)*exp(-b[1]-(x%*%b[-1])))+lambda*sum(b[-1]^2));
			} 
		  if (surrogate == 'sqhinge'){  
			return (sum(y1*(y1>0)*(pmax(1-b[1]-(x%*%b[-1]),0))^2)+
			sum(-y1*(y1<0)*(pmax(1+b[1]+(x%*%b[-1]),0))^2)+
			sum(y2*(y2>0)*(pmax(1+b[1]+(x%*%b[-1]),0))^2)+
			sum(-y2*(y2<0)*(pmax(1-b[1]-(x%*%b[-1]),0))^2)+lambda*sum(b[-1]^2));
			} 	 		
	      if (surrogate == 'hinge'){  
			return (sum(y1*(y1>0)*(pmax(1-b[1]-(x%*%b[-1]),0)))+
			sum(-y1*(y1<0)*(pmax(1+b[1]+(x%*%b[-1]),0)))+
			sum(y2*(y2>0)*(pmax(1+b[1]+(x%*%b[-1]),0)))+
			sum(-y2*(y2<0)*(pmax(1-b[1]-(x%*%b[-1]),0)))+lambda*sum(b[-1]^2));
			} 	 			
		}
    #bTilde = optim (rnorm(p+1), objFn, method="BFGS")$par;
    bTilde = optim (rep(1,p+1), objFn, method="BFGS")$par;
	#bTilde = bTilde/abs(bTilde[1]);
	return(list(bTilde = bTilde, wpos1 = y1, wneg1 = y2))
   
}
 