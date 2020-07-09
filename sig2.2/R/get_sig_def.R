# The function to get the failure time of a failure permutation
# x: 1 possible failure permutation
fail_num<-function(x,cutsets)
{
    for(i in 1:length(x))
    {
        for(cutset in cutsets)
        {
            # if the first i components in the failure permutation contains a cutset, then break and return i
            if(sum(x[1:i] %in% cutset)==length(cutset)) return(i)
        }
    }
}

# Use Da`s method to compute the signature and survival signature of an arbitrary system
# Define a graph of the system using the igraph package, with additional nodes named ``s'' and ``t'' for either end of the system diagram
get_sig_def<-function(graph=NULL,cutsets=NULL,frac=TRUE)
{
    ptm <- proc.time()
    if(length(cutsets)==0) {
        if(sum(is.na(match(c("s","t"), igraph::V(graph)$name)))) {
            stop("Oops! It seems that you didn't give us anything or you gave us something in a wrong way! We can't help you compute the signature!")}

        cutsets=minimalVertexPairCutSets(graph,'s','t')
        # unlist(cutsets[1]),you will get some numbers named by elements' names
        # give these numbers to igraph::V(graph)$name[] you will get the elements' names
        # so the cutsets(used to be a graph list) now turn to be a list including elements' names
        cutsets=lapply(cutsets, function(x) { igraph::V(graph)$name[unlist(x)] })
        # n is the number of whole graph(including s and t)
        n=length(igraph::V(graph)$name)
    }

    # get all components in this system and the number of components
    comp=c()
    for(i in 1:length(cutsets)){
        comp=union(comp,unlist(cutsets[i]))
    }
    totalC=length(comp)
    n=totalC+2

    permn_sig=tabulate(unlist(combinat::permn(comp,fail_num,cutsets)),nbins = n-2)
    #given that the 1st argument of tabulate() must be a a numeric vector (of positive integers),
    #or a factor,while the result of combinat::permn() is a list, so we use the unlist to get a simple vector
    #including the failure times of each permutation
    sig<-permn_sig/factorial(n-2)

    # get survSig from sig
    survSig<-c()
    for(i in 1:(length(sig)-1)){
        survSig<-append(survSig,sum(sig[(i+1):length(sig)]))
    }
    survSig=append(survSig,0)

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

