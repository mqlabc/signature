# The function to get the great greatest common factor of two nummbers
single_gcd<-function(a,b)
{
    if(a==0) return(b)
    else
    {
        if(a<b)
        {
            t=b;b=a;a=t
        }
        q<-b
        a<-a%%b
        single_gcd(a,q)
    }
}

# The function to get the great greatest common factors of corresponding 2 components in 2 vectors
gcd<-function(a,b)
{
    c=c()
    if(length(a)!=length(b)) print('Wrong! When getting the gcds, 2 vectors have different lengths!')
    else
    {
        for(i in 1:length(a))
            c=append(c,single_gcd(a[i],b[i]))
    }
    return(c)
}

# is.wholenumber
is.wholenumber<-function (x, tol = .Machine$double.eps^0.5)
    abs(x - round(x)) < tol

# the fractionizing function to be called
fra<-function (x, j = 7)
{
    i = 0
    y = seq(j)
    a = seq(j)
    for (q in 1:j) {
        y[q] = x * (10^q)
    }
    if(x<1e-10){
        return('0')
    }
    if(abs(x-1)<1e-10){
        return('1')
    }
    if (is.wholenumber(y[j]) == TRUE) {
        k = round(y[j])
        b = gcd(k, (10^(j)))
        k = k/b
        b = (10^(j))/b
        c = paste(k, "/", b)
        return(c)
        i = 1
    }
    for (q in 1:j) {
        a[q] = y[q] - x
    }
    a[7] = round(y[7])
    if (is.wholenumber(a[1]) == TRUE & j > 1) {
        a[1] = round(a[1])
        b = gcd(a[1], 9 * 10^(7 - j))
        a[1] = a[1]/b
        b = 9 * 10^(7 - j)/b
        c = paste(a[1], "/", b)
        return(c)
        i = 1
    }
    else if (is.wholenumber(a[2]) == TRUE & j > 2) {
        a[2] = round(a[2])
        b = gcd(a[2], 99 * 10^(7 - j))
        a[2] = a[2]/b
        b = 99 * 10^(7 - j)/b
        c = paste(a[2], "/", b)
        return(c)
        i = 1
    }
    else if (is.wholenumber(a[3]) == TRUE & j > 3) {
        a[2] = round(a[3])
        b = gcd(a[3], 999 * 10^(7 - j))
        a[3] = a[3]/b
        b = 999 * 10^(7 - j)/b
        c = paste(a[3], "/", b)
        return(c)
        i = 1
    }
    else if (is.wholenumber(a[4]) == TRUE & j > 4) {
        a[4] = round(a[4])
        b = gcd(a[4], 9999 * 10^(j - 7))
        a[4] = a[4]/b
        b = 9999 * 10^(7 - j)/b
        c = paste(a[4], "/", b)
        return(c)
        i = 1
    }
    else if (is.wholenumber(a[5]) == TRUE & j > 5) {
        a[5] = round(a[5])
        b = gcd(a[5], 99999 * 10^(7 - j))
        a[5] = a[5]/b
        b = 99999 * 10^(7 - j)/b
        c = paste(a[5], "/", b)
        return(c)
        i = 1
    }
    else if (is.wholenumber(a[6]) == TRUE & j > 6) {
        a[6] = round(a[6])
        b = gcd(a[6], 999999 * 10^(7 - j))
        a[6] = a[6]/b
        b = 999999 * 10^(7 - j)/b
        c = paste(a[6], "/", b)
        return(c)
        i = 1
    }
    if (i == 0 & j > 1) {
        j = j - 1
        x = 10 * x
        return(fra(x, j))
    }
    if (i == 0 & j == 1) {
        k = round(10 * x)
        b = gcd(k, 1e+08)
        k = k/b
        b = 1e+08/b
        c = paste(k, "/", b)
        j = 0
        return(c)
    }
}

