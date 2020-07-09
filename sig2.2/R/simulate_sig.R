# unrecommended function call because of low speed
get_num_from_graph<-function(graph,m){
    n=length(igraph::V(graph)$name)-2
    # fail_num
    fail_num=c()

    for(k in 1:m)
    {
        # Use the runif() function to generate random numbers on the interval (0,1) whose quantity is identical to the number of elements
        # Use rank() to get the order from small to large which is 1,2,3...(rank statistics)
        # fail_perm() is the failure order of simulated components
        fail_perm=rank(stats::runif(n))
        j=0
        graph1=graph

        # Calculate the number of failure components which cause the system to fail according the corresponding fail_perm when the graph is connected
        # This is the system death time
        while(igraph::edge_connectivity(graph1)!=0)
        {
            j=j+1
            # Remove the nodes from the original image
            graph1=igraph::delete.vertices(graph,as.character(fail_perm[1:j]))
        }
        # Each k is a simulation, and the resulte j is what we care about
        # fail_num is a vector with these values as components
        fail_num=append(fail_num,j)
    }
    # Calculate the system death time of each simulation

    return(fail_num)
}

get_num_from_cutsets<-function(cutsets,m){
    #each component in fail_num stands for at the jth attack, the system is dead
    fail_num=c()
    u=c()
    for (i in 1:length(cutsets)) {
        u=union(u,unlist(cutsets[i]))
    }
    n=length(u)

    for(k in 1:m){
        fail_perm=rank(stats::runif(n))
        j=1
        sign=0
        # sign equals to zero means that the system is not fail
        while(sign==0){
            for(cutset in cutsets){
                if(sum(u[fail_perm[1:j]]%in%cutset)==length(cutset)){
                    sign=1
                    break
                }
            }
            j=j+1
        }
        fail_num=c(fail_num,j-1)
    }
    return(list(fail_num,n))
}

# Signature is calculated by simulation method(defaulted 10000 times)
simulate_sig<-function(graph=NULL,cutsets=NULL,m=10000)
{
    ptm <- proc.time()
    if(length(cutsets)==0){
        if(sum(is.na(match(c("s","t"), igraph::V(graph)$name)))) {
            stop("Oops! It seems that you didn't give us anything or you gave us something in a wrong way! We can't help you compute the signature!")}

        cutsets=minimalVertexPairCutSets(graph,'s','t')
        cutsets=lapply(cutsets, function(x) { igraph::V(graph)$name[unlist(x)] })

        totalC=length(igraph::V(graph)$name)-2
        fail_num=unlist(get_num_from_cutsets(cutsets,m)[1])
        #fail_num=get_num_from_graph(graph,m)
    }else{
        result=get_num_from_cutsets(cutsets,m)
        totalC=unlist(result[2])
        fail_num=unlist(result[1])
    }
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

    results=data.frame(No.=1:totalC,EstimatedSignature=sig,EstimatedSurvivalSignature=survSig,Time_used=(proc.time()-ptm)[[3]])
    return(results)
}

