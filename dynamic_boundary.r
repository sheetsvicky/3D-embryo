# 2015.09.29, calculate the clonal boudary from 26-cell stage to the end of the gastrulation. Recover the script from "dynamic_boundary.rdata" of last week from a lost because of airdriver problem.

founder<-c('ABala','ABalp','ABara','ABarp','ABpla','ABplp','ABpra','ABprp','MS','E','C','D','P') # 13
founder_index<-1:13
names(founder_index)<-founder
founder_col<-c('navy','plum','cyan','magenta','darkgreen','khaki','chartreuse','yellow','grey','sienna','tomato','lightseagreen','black','white')
names(founder_col)<-c(founder,'ud')

# calculate space distance, x y z ==> 6-8 columns
trans_z<-0.254 ### !!!Z is the plane number, and x, y, diameter are in pixels, all Zhuo's image are the same, check the XML file for this
cal_euclidean<-function(x,data){
  return(sqrt(apply((matrix(rep(as.numeric(x),dim(data)[1]),nc=3,byrow=T)-data[,6:8])^2,1,sum)))
}

# get the ancestor of nuclei by name
get_ancestor<-function(x){
  res<-'ud' # undetermined
  if(nchar(x)==0|substring(x,0,1)=="N") return(res)
  else if(substring(x,0,1)!="A"){
    res<-substring(x,0,1)
    if(res=='M') res<-'MS'
    if(res=='Z') res<-'P'
  }else if(nchar(x)>4){
    res<-substring(x,0,5)
  }else res<-x
  return(res)
}

# calculate the dynamic at each time point
dynamic_ashape<-
function(emb,end_time,alpha_value=0.4,is_cell=T,num=1:13){ # from 26-cell to 350-cell time determined by lineage tree
  emb_path<-"./embryos_data/"
  i<-1
  while(T){ # determine the time for 26-cell stage
    time_txt<-i
    if(i<10) time_txt<-paste('0',time_txt,sep='')
    if(i<100) time_txt<-paste('0',time_txt,sep='')
    nu_tm_raw<-read.table(paste(emb_path,emb,'/nuclei/t',time_txt,'-nuclei',sep=''),sep=',',as.is=T,strip.white=T) 
    nu_tm_valid<-nu_tm_raw[nu_tm_raw[,2]==1&sapply(nu_tm_raw[,10],substring,0,1)!='N',]
    if(dim(nu_tm_valid)[1]>=26) break
    i<-i+1
  }
  start_time<-i
  time_res<-lapply(start_time:end_time,function(time){
    time_txt<-time
    if(time<10) time_txt<-paste('0',time_txt,sep='')
    if(time<100) time_txt<-paste('0',time_txt,sep='')
    nu_tm_raw<-read.table(paste(emb_path,emb,'/nuclei/t',time_txt,'-nuclei',sep=''),sep=',',as.is=T,strip.white=T) 
    nu_tm_valid<-nu_tm_raw[nu_tm_raw[,2]==1&sapply(nu_tm_raw[,10],substring,0,1)!='N',]
    nu<-nu_tm_valid[,10]
    nu_an<-sapply(nu,get_ancestor)
    nu_an_order<-nu_an[order(founder_index[nu_an],nu)] # order according to their ancestor and their AP position
    nu_tm_order<-nu_tm_valid[order(founder_index[nu_an],nu),]
    nu_tm_order[,8]<-nu_tm_order[,8]/trans_z 
    nu_dist<-apply(nu_tm_order[,6:8],1,cal_euclidean,data=nu_tm_order)
    ap<-max(nu_dist)
    nu_dist<-nu_dist/max(nu_dist) # normalize
#    rownames(nu_dist)<-nu_tm_order[,10]->colnames(nu_dist)
    if(is_cell){
      ashape_cell<-lapply(founder[num],function(x){
        if(sum(nu_an_order==x)<=4){ # do not model ashape when cell number is less than 4
          bound_cell<-as.matrix(cbind(nu_tm_order[nu_an_order==x,6:8],rep(1,sum(nu_an_order==x)),matrix(rep(2,length(alpha_value)*sum(nu_an_order==x)),nc=length(alpha_value))))
        }else{
	  points<-as.matrix(nu_tm_order[nu_an_order==x,6:8])
	  if(length(unique(points[,3]))==1) points[1,3]<-points[1,3]+0.1 # if all cells are coplanar, add a bit error on Z plane
          ashape<-ashape3d(points,alpha=ap*alpha_value)
          bound_cell<-cbind(ashape$x,ashape$vertex[,c(-1,-3,-4)]) # return the XYZ, whether the cell belong to convex hull, and how the cell belongs to alpha-complex (which is slightly different to alpha-shape, 1,2 or 3 indicates interior, regular or singular)
        }
boundary_index<-rownames(bound_cell)[bound_cell[,5]==2]
    lineage<-unique(nu_an[rownames(nu_tm_valid) %in% rownames(bound_cell)]) # this step is how to get the ancestor
    if(length(lineage)!=1) print("Lineage not correct!")
    clonal_distance<-nu_dist[nu_an_order==lineage,nu_an_order!=lineage]
if(sum(bound_cell[,5]==2)==1){ # deal with only one cell
  bound_distance<-matrix(clonal_distance,nr=1) 
  colnames(bound_distance)<-names(clonal_distance)
  rownames(bound_distance)<-rownames(bound_cell)
}else{ 
  bound_distance<-clonal_distance[rownames(clonal_distance) %in% boundary_index,]    
}
    min_cell<-sapply(nu_tm_valid[colnames(bound_distance)[apply(bound_distance,1,which.min)],10],get_ancestor)
names(min_cell)<-rownames(bound_distance)
min_cell_res<-rep(0,dim(bound_cell)[1])
names(min_cell_res)<-rownames(bound_cell)
min_cell_res[names(min_cell)]<-sapply(min_cell,function(y){which(founder==y)})
bound_cell_res<-cbind(bound_cell,min_cell_res)
colnames(bound_cell_res)<-c('X','Y','Z','chull','ashape','closest lineage')
return(bound_cell_res)
      })
      return(ashape_cell)
    }else{
      ashape<-lapply(founder[num],function(x){
        if(sum(nu_an_order==x)==4){ # When cell number is equal to 4, model as a convex hull
   points<-as.matrix(nu_tm_order[nu_an_order==x,6:8])
   if(length(unique(points[,3]))==1) points[1,3]<-points[1,3]+0.1 # if all cells are coplanar, add a bit error on Z plane
   return(list(nu_tm_order[nu_an_order==x,6:8],convhulln(points))) 
}else if(sum(nu_an_order==x)<4){ # When cell number is less than 4, do not model, just return the vertex
   return(list(nu_tm_order[nu_an_order==x,6:8]))
}else{
	points<-as.matrix(nu_tm_order[nu_an_order==x,6:8])
	if(length(unique(points[,3]))==1) points[1,3]<-points[1,3]+0.1 # if all cells are coplanar, add a bit error on Z plane
        return(ashape3d(points,alpha=ap*alpha_value))
}
      })
      return(ashape) # return the whole ashape object
    }
  })
  names(time_res)<-start_time:end_time
  return(time_res)
}

# plot embryo in lineage
view_emb_ashape<-
function(emb,clone,which_alpha=1,open_new=T,is_print=F,is_save=F,dir='plot/',time='0',view_mx=par3d('userMatrix')){
  if(open_new){
     if(is_save){
       par3d(windowRect=c(0,0,700,700)) 
#    par3d(viewport=c(0,0,1200,1200),windowRect=c(0,0,600,600))
     }else{
       open3d()
       bg3d("white")
     }
  }else{
    if(is_save) par3d(windowRect=c(0,0,700,700)) 
    rgl.clear()
  }  
#  title3d(main = paste("TTT")) # cannot work
  sapply(clone,function(x){
    if(length(emb[[x]])>2){ 
       plot(emb[[x]],col=rep(founder_col[x],3),clear=F,transparency=0.8,indexAlpha=which_alpha)
       if(is_print) cat(founder[x],"a-shape","\n")
    }else if(length(emb[[x]])==2){ # plot convex hull
       index<-as.vector(t(emb[[x]][[2]]))
       rgl.triangles(emb[[x]][[1]][,1][index],emb[[x]][[1]][,2][index],emb[[x]][[1]][,3][index],col = founder_col[x],transparency=0.8)
       if(is_print) cat(founder[x],"chull","\n")
    }else if(dim(emb[[x]][[1]])[1]==3){ # three vertex, plot a triangle
       rgl.triangles(emb[[x]][[1]][,1],emb[[x]][[1]][,2],emb[[x]][[1]][,3],col = founder_col[x],transparency=0.8)
       if(is_print) cat(founder[x],"triangle","\n")       
    }else if(dim(emb[[x]][[1]])[1]==2){ # two vertex, plot a line
       rgl.lines(emb[[x]][[1]][,1],emb[[x]][[1]][,2],emb[[x]][[1]][,3],col = founder_col[x],transparency=0.8,lwd=2)
       rgl.points(emb[[x]][[1]][,1],emb[[x]][[1]][,2],emb[[x]][[1]][,3],col = founder_col[x],transparency=0.8,size=10)
       if(is_print) cat(founder[x],"line","\n")       
    }else{ # one vertex, plot a point
       rgl.points(emb[[x]][[1]][,1],emb[[x]][[1]][,2],emb[[x]][[1]][,3],col = founder_col[x],transparency=0.8,size=10)
       if(is_print) cat(founder[x],"point","\n")       
    }
  })
  if(is_save){
     rgl.viewpoint(userMatrix=view_mx) # adjust the view
     time_txt<-time
     if(as.numeric(time)<100) time_txt<-paste('0',time_txt,sep='')
     file_name<-paste(dir,'t',time_txt,'.png',sep='')
     rgl.snapshot(filename=file_name,fmt='png')
     rgl.close()
  }
  return(par3d("userMatrix"))
}
 
# calculate the first WT sample
emb_info<-read.table('./embryo_lineage_refine_20150812.csv',sep=',',fill=T,comment.char='',head=T,as.is=T,row.name=1) # 207  6, contains the development stage, i.e., the end time
wt_sample<-c('emb1','emb2')
wt_ashape<-dynamic_ashape(emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],is_cell=F) # a list of list of ashape 
wt_ashape_cell<-dynamic_ashape(emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],is_cell=T) # a list of list of matrices
# view in 3D
view_emb_ashape(wt_ashape[[12]],clone=1:13)
# view all time point
sapply(wt_ashape,view_emb_ashape,clone=1:13,is_save=T,dir='wt_sample1')
# To gif by ImageMagick, install by macports
convert -delay 30 *.png -loop 0 wt_sample1.gif


# calculate connectedness based probability and the closest cell
cal_prob_connect<-
function(ashape_cell,num=12){ # for one time point
  connect<-sapply(ashape_cell,function(x){
    min_vec<-rep(0,13)
    names(min_vec)<-1:13
    min_num<-tapply(x[x[,5]==2,6],x[x[,5]==2,6],length)
    min_vec[names(min_num)]<-min_num
    return(min_vec)
  })
  lin_num<-sapply(ashape_cell,function(x){dim(x)[1]})
  bou_num<-sapply(ashape_cell,function(x){sum(x[,5]==2)})
  total<-sum(lin_num)
  is_con<-matrix(rep(0,num*num),nr=num) # a matrix, whose lower triangle store whether two clonal boundary are connected.
  rownames(is_con)<-founder[1:num]->colnames(is_con)
  for(i in 1:num){
    for(j in 1:i){
      if(connect[i,j]>=(lin_num[i]*bou_num[j]/(total-lin_num[j]))|connect[j,i]>=(bou_num[i]*lin_num[j]/(total-lin_num[i]))) is_con[i,j]<-1
    }
  }
  return(is_con)
}

wt_connect<-lapply(wt_ashape_cell,cal_prob_connect)

# compare connected
compare_connect<-
function(x,y,num=12,is_print=F){ 
  total<-choose(num,2)
  diff<-0
  for(i in 1:num){
    for(j in 1:i){
      if((x[i,j]+y[i,j])==1){
         if(is_print) cat(rownames(x)[i],"-",rownames(x)[j],"\t",x[i,j],"\t",y[i,j],"\n",sep='')
 diff<-diff+1
      }
    }
  }
  return(diff*100/total)
}

# print connectedness
print_connect<-
function(x,num=12){ # default: ignore P 
  for(i in 1:num){ 
    cat(rownames(x)[i],": ",sep='')
    for(j in 1:num){
      if(i<j){
        index1<-j
index2<-i
      }else{
        index1<-i
index2<-j
      }
      if(x[index1,index2]==1) cat(rownames(x)[j],",",sep='')
    }
    cat("\n")
  }
}

# The cell number in each time point
sapply(wt_ashape_cell,function(x){sum(sapply(x,function(y){dim(y)[1]}))})
# The difference of connectedness compared to last time point
plot(1:120,sapply(wt_connect,compare_connect,y=wt_connect[[120]])) # 10% difference with 280-cell stage
# The number of connectedness in each time point
plot(1:120,sapply(wt_connect,sum))

