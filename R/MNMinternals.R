Ginv<-function(A)
{
p<-sqrt(dim(A)[1])
eig <- eigen(A)
ind <- 1:((p+2)*(p-1)/2)
eig$values[ind] <- 1/eig$values[ind]
res<-eig$vectors %*% tcrossprod(diag(eig$values), eig$vectors)
Re(res)
}

Cpp<-function(p)
{
I <- diag(p)
J <- tcrossprod(as.vector(I))
K<-matrix(0,ncol=p^2,nrow=p^2)
    for (i in 1:p)
        {
        for (j in 1:p)
            {
             K <- K + kronecker(outer(I[,i],I[,j]), outer(I[,j],I[,i]))
            }
        }
(diag(p^2)+K)/2-J/p
}

RowNorms <- function(X)
    {
    sqrt(rowSums(X^2))
    }
