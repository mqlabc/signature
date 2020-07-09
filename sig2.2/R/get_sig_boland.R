# convert x to the binary format
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
# Define a graph of the system using the igraph package, with additional nodes named ``s'' and ``t'' for either end of the system diagram
get_sig_boland <- function(graph=NULL, cutsets = NULL, frac = TRUE) {
    ptm <- proc.time()
    if(length(cutsets)==0) {
        if(sum(is.na(match(c("s","t"), igraph::V(graph)$name)))) {
            stop("Oops! It seems that you didn't give us anything or you gave us something in a wrong way! We can't help you compute the signature!")}

        cutsets <- minimalVertexPairCutSets(graph, "s", "t")
        cutsets=lapply(cutsets, function(x) { igraph::V(graph)$name[unlist(x)] })
    }

    # get all components in this system and the number of components
    comp=c()
    for(i in 1:length(cutsets)){
        comp=union(comp,unlist(cutsets[i]))
    }
    totalC=length(comp)

    # get the number of cutsets of size1
    size1_cut=0
    for(cutset in cutsets){
        if(length(cutset)==1){
            size1_cut=size1_cut+1
        }
    }


    # Now start finding the survival signature
    # size of each path set
    size_of_eps=c()
    for(x in 0:{2^totalC-1}) {
        #!!!!! change x to binary????
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
    # r(n-i) is equal to raw_numerator; the 1st element of rev_tab corr espond the survSig[0], which means nothing
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
    survSig1<-survSig[2:length(survSig)]
    sig_from_2<-survSig[1:(length(survSig)-1)]-survSig1
    # add the first signature component
    sig<-c(size1_cut/totalC,sig_from_2)

    #如果是小数
    if(!frac){
        results=data.frame(No.=1:totalC,Signature=sig,SurvivalSignature=survSig,Time_used=(proc.time()-ptm)[[3]])
    }else{
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
        results=data.frame(No.=1:totalC,Signature=sig_str,SurvivalSignature=survSig_str,Time_used=(proc.time()-ptm)[[3]],stringsAsFactors = FALSE)
    }
    return(results)
}
