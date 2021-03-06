# Convert x to the binary format
digitsBase<-function (x, base = 2, ndigits = 1 + floor(1e-09 + log(max(x,1), base)))
{
    if (any(x < 0))
        stop("'x' must be non-negative integers")
    if (any(x != trunc(x)))
        stop("'x' must be integer-valued")
    r <- matrix(0, nrow = ndigits, ncol = length(x))
    if (ndigits >= 1)
        for (i in ndigits:1) {
            r[i, ] <- x%%base
            if (i > 1)
                x <- x%/%base
        }
    class(r) <- "basedInt"
    attr(r, "base") <- base
    r
}


# Use Boland`s method to compute the signature and survival signature of an arbitrary system
# Define a graph of the system using the igraph package, with additional nodes named "s" and "t" for either end of the system diagram
get_survSig_boland <- function(totalC, cutsets, num_of_cutsets) {
    # Now start finding the survival signature
    # to get size of each path set
    size_of_eps=c()
    for(x in 0:{2^totalC-1}) {
        # change x to binary
        state <- digitsBase(x, base=2, ndigits=totalC)
        failed <- comp[state==0]
        # True is a necessary argument in vapply, 遍历到各节点,但以cutset为单位判断某cutset的所有元素是否都出现再failed(都是编号而非节点名)里,如果都在(prod==1)该cutset返回True,全部cutset对应都返回False则系统工作
        if( sum(vapply(cutsets, function(x) { prod(x %in% failed) == 1 }, TRUE)) == 0 ) { # TRUE only if the system is working
            size_of_ps=sum(state)
            size_of_eps<-append(size_of_eps,size_of_ps)
        }
    }
    tab=tabulate(size_of_eps)
    # rev_tab is the reverse version of tab
    rev_tab=rev(tab)
    # generrate the r(n-i), which means the nunmber of the path sets of size n-i
    # r(n-i) is equal to raw_numerator; the 1st element of rev_tab correspond the survSig[0], which means nothing
    raw_numerator=c()
    for(i in 2:length(rev_tab)){
        raw_numerator=append(raw_numerator,rev_tab[i])
    }

    # the denominator
    raw_denominator=c()
    for(i in 1:(totalC-1)){
        raw_denominator=append(raw_denominator,choose(totalC,i))
    }

    survSig<-c(raw_numerator/raw_denominator,0)
    survSig
}










# Get the size of the union of cutsets
count_union_num<-function(cutsets,colu,totalC){
    vec=c()
    ele=cutsets[colu]

    for(i in ele){
        vec=union(as.vector(unlist(i)),vec)
    }
    # since the last component of survival signature is 0, it`s useless to know the union of size totalC
    return(length(vec))
}

# list[1]: get the number of components in each possible union of different cutsets
# list[2]: the corresponding number of cutsets
get_list<-function(cutsets,totalC){
    num=c()
    symbol=c()

    # The algebraic union of all minimal cut sets is the set of all the system’s components
    # obviously, the last component of survival signature is 0
    # so length(cutsets)-1 can filter some cases which causes the union of cutsets contains all components
    for(i in 1:(length(cutsets)-1)){
        # combn() returns a matrix where a column represents a combnation of size i with the components in 1:length(cutsets)
        # 1:?
        comb_mat=combinat::combn(length(cutsets),i)
        for(j in 1:ncol(comb_mat)){
            # comb_mat[,j] is the indexes of the selected cutsets
            num=append(num,count_union_num(cutsets,comb_mat[,j],totalC))
            # i is the corresponding number of cutsets
            symbol=append(symbol,i)
        }
    }
    return(list(num,symbol))
}


# l2 is the 2nd element of get_list(); return 1 if the number of cutsets is even, else -1
to_symbol<-function(l2){
    v=c()
    for(i in l2){
        if(i%%2==0){
            n=1
        }else{
            n=-1
        }
        v=append(v,n)
    }
    return(v)
}

# called by get_survSig_m(); totalC is the number of all components
# return the 'left multiply factor' consisting of ratios; m is the index of survSig with a range of [1,totalC-1]
get_lmf<-function(l1,totalC,m){
    # l1<=m guarantees selecting the union size <=m
    nd=l1[l1<=m]
    # length of numerator depends on nd
    numerator=choose(m,nd)
    # length of denominator depends on nd
    denominator=choose(totalC,nd)

    return(numerator/denominator)
}

# input the list, totalC and the mth index of survial signature
# return the mth component
get_survSig_m<-function(l,totalC,m){
    lmf=get_lmf(unlist(l[1]),totalC,m)
    # converted to 1 or -1 according to the number of cutsets union
    #important!!!!
    symb=to_symbol(unlist(l[2]))[unlist(l[1])<=m]
    ssig_m=1+as.numeric(lmf%*%symb)
    return(ssig_m)
}

# Use Da`s method to compute the survival signature of an arbitrary system
get_survSig_da<-function(totalC, cutsets, num_of_cutsets){
    # initialise survSig to be a vector of zeros?????????????????????????????? ones() or what?
    survSig<-rep(0,totalC)
    # If not a parallel system
    if(num_of_cutsets!=1){
        l=get_list(cutsets,totalC)

        # get the [1:totalC-1] components ofsurvSig
        for(i in 1:(totalC-1)){
            ss=get_survSig_m(l,totalC,i)
            if(ss<1e-8){
                # if there is a 0 in survSig, the components behind it are all 0s
                break
            }else{
                survSig[i]<-ss
            }
        }
    }else{
        survSig[1:(totalC-1)]=1
    }
    survSig
}










get_num_from_cutsets<-function(cutsets,comp,m){
    #each component in fail_num stands for at the jth attack, the system is dead
    fail_num=c()
    totalC=length(comp)

    for(k in 1:m){
        fail_perm=rank(stats::runif(totalC))
        j=1
        sign=0
        # sign equals to zero means that the system is not fail
        while(sign==0){
            for(cutset in cutsets){
                if(sum(comp[fail_perm[1:j]]%in%cutset)==length(cutset)){
                    sign=1
                    break
                }
            }
            j=j+1
        }
        fail_num=c(fail_num,j-1)
    }
    fail_num
}

# Signature is calculated by simulation method(defaulted 10000 times)
simulate_sig<-function(totalC, cutsets, comp, m=10000)
{
    fail_num=get_num_from_cutsets(cutsets,comp = comp,m)
    stat=tabulate(fail_num)


    # After the system fail, each component becomes a zero, and the number of zero components is n_0, according to the definition of signature
    n_0=totalC-max(fail_num)
    # The zero will not appear before the non-zero component of signature
    sig=append(stat/m,rep(0,times=n_0))

    # get survSig from sig
    survSig<-c()
    for(i in 1:(length(sig)-1)){
        survSig<-append(survSig,sum(sig[(i+1):length(sig)]))
    }
    survSig=append(survSig,0)

    results=data.frame(No.=1:totalC,EstimatedSignature=sig,EstimatedSurvivalSignature=survSig)
    results
}












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

output_form<- function(survSig,size1_cut,totalC,time_used,frac=TRUE) {
    survSig1<-survSig[2:length(survSig)]
    sig_from_2<-survSig[1:(length(survSig)-1)]-survSig1
    # Add the first signature component
    sig<-c(size1_cut/totalC,sig_from_2)
    # survSig is of length toatalC
    totalC=length(survSig)
    # If wanted in decimal
    if(!frac){
        results=data.frame(No.=1:totalC,Signature=sig,SurvivalSignature=survSig,Time_used=time_used)
    }else{
        # if wanted in decimal
        sig_str=c()
        for(i in sig){
            s=gsub(' ','',fra(i))
            sig_str=c(sig_str,s)
        }

        survSig_str=c()
        for(i in survSig){
            s=gsub(' ','',fra(i))
            survSig_str=c(survSig_str,s)
        }
        results=data.frame(No.=1:totalC,Signature=sig_str,SurvivalSignature=survSig_str,Time_used=time_used,stringsAsFactors = FALSE)
    }
    results
}









# Find collection of minimal (s,t) vertex cut sets
minimalVertexPairCutSets <- function(graph, start, terminate) {
    # Get all separators
    stCutSets <- igraph::minimal.st.separators(graph)
    # Double check separator (there was a bug in igraph early 0.6 releases)
    stCutSets <- subset(stCutSets, vapply(stCutSets, igraph::is.minimal.separator, TRUE, graph=graph))
    # Obviously, if s or t are in the cut set is it not of interest
    stCutSets <- stCutSets[vapply(lapply(stCutSets, match, match(c(start, terminate), igraph::V(graph)$name), 0), sum, 1)==0]
    # Check that the cut set *is* an s,t separator, since igraph gets all separators
    cutGraphs <- lapply(stCutSets, igraph::delete.vertices, graph=graph)
    reachable <- lapply(cutGraphs, igraph::subcomponent, start, "out")
    isSTsep <- rep(TRUE, length(stCutSets))
    for(j in 1:length(stCutSets)) { # Can't do nice vectorised stuff, because when deleting vertices igraph changes all the ids, so each check different
        if(match(terminate, igraph::V(cutGraphs[[j]])$name) %in% reachable[[j]]) {
            isSTsep[j] <- FALSE
        }
    }
    if(sum(isSTsep) < 1) return(NULL);
    return(stCutSets[isSTsep])

}
auto_judge<-function(graph=NULL,cutsets=NULL,frac=TRUE){
    ptm <- proc.time()
    if(length(cutsets)==0) {
        if(sum(is.na(match(c("s","t"), igraph::V(graph)$name)))) {
            stop("Oops! It seems that you didn't give us anything or you gave us something in a wrong way! We can't help you compute the signature!")}

        cutsets=minimalVertexPairCutSets(graph,'s','t')
        # converted to the form of list, each component is a cutset consisting of str of comp's name
        cutsets=lapply(cutsets, function(x) { igraph::V(graph)$name[unlist(x)] })
        # n is the number of whole graph( including s and t)
        n=length(igraph::V(graph)$name)
    }

    # get the n==totalC+2
    # comp: all unique components in this system and the number of components
    comp=c()
    for(i in 1:length(cutsets)){
        comp=union(comp,unlist(cutsets[i]))
    }
    totalC=length(comp)

    # get the number of cutsets and the number of size1_cutsets
    num_of_cutsets=length(cutsets)
    size1_cut=length(cutsets[lengths(cutsets)==1])

    if(totalC>30 | num_of_cutsets>30){
        # simulate
        no_time=simulate_sig(totalC, cutsets, comp)
        results=data.frame(no_time,Time_used=(proc.time()-ptm)[[3]])
    }else if(totalC >= num_of_cutsets){
        # da's
        results=output_form(get_survSig_da(totalC, cutsets, num_of_cutsets), size1_cut, totalC, time_used = (proc.time()-ptm)[[3]], frac = TRUE)
    }else{
        # boland
        results=output_form(get_survSig_boland(totalC, cutsets, num_of_cutsets), size1_cut, totalC, time_used = (proc.time()-ptm)[[3]], frac = TRUE)
    }
}
