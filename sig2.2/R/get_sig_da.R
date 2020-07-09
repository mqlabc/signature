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

get_survSig<-function(graph=NULL,cutsets=NULL){
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
    n=totalC+2

    # get the number of cutsets and the number of size1_cutsets
    num_of_cutsets=length(cutsets)
    size1_cut=0
    for(cutset in cutsets){
        if(length(cutset)==1){
            size1_cut=size1_cut+1
        }
    }

    # initialise survSig to be a vector of zeros?????????????????????????????? ones() or what?
    survSig<-stats::runif(totalC,min = 0,max = 0)
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
    return(list(survSig,totalC,size1_cut))
    # survSig is of length toatalC
}

# Use Da`s method to compute the signature and survival signature of an arbitrary system
# Define a graph of the system using the igraph package, with additional nodes named ``s'' and ``t'' for either end of the system diagram
get_sig_da<-function(graph=NULL,cutsets=NULL,frac=TRUE){
    ptm <- proc.time()
    li=get_survSig(graph,cutsets)
    survSig<-unlist(li[1])
    totalC<-unlist(li[2])
    size1_cut<-unlist(li[3])

    survSig1<-survSig[2:length(survSig)]
    sig_from_2<-survSig[1:(length(survSig)-1)]-survSig1
    # add the first signature component
    sig<-c(size1_cut/totalC,sig_from_2)


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
