# Functions common to most models
##################################

# comp.gt
#########
# logical function to see if all elements of x are > u
# used via "apply" on matrices of data to determine which rows correspont to completely uncensored observations

comp.gt<-function(x,u){all(x>u)}

# BCi
#####

# Inverse Box--Cox transformation, using limit form if gamma close to 0

BCi<-function(x,gamma,sig)
{
  y<-NULL
  d<-length(x)
  for(j in 1:d)
  {  
    if(abs(gamma[j])>1e-6)
    {
      y[j]<-(1+gamma[j]*x[j]/sig[j])^(1/gamma[j])
    }
    else{
      y[j]<-exp(x[j]/sig[j])
    }
  }
  return(y)
}

# Jac
#####

# Jacobian for transformation from multivartiate Pareto to MGP, using limit for gamma if near 0

Jac<-function(x,gamma,sig)
{
  J<-NULL
  d<-length(x)
  for(j in 1:d)
  {  
    if(abs(gamma[j])>1e-6)
    {
      J[j]<-(1/sig[j])*(1+gamma[j]*x[j]/sig[j])^(1/gamma[j]-1)
    }
    else{
      J[j]<-(1/sig[j])*exp(x[j]/sig[j])
    }
  }
  return(prod(J))
}
