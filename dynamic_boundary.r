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
  emb_path<-"/Volumes/baolab/xuy/embryos_data/"
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
emb_info<-read.table('/Volumes/baolab/xuy/distance_mx/lineage_loss/embryo_lineage_refine_20150812.csv',sep=',',fill=T,comment.char='',head=T,as.is=T,row.name=1) # 207  6, contains the development stage, i.e., the end time
wt_sample<-c("ZD_RW10348_WT_20110126_2_s2_emb1_edited", "ZD_BV82_WT_20110426_1_s1_emb3_edited","ZD_RW10714_WT_20110724_1_s3_emb3_edited", "ZD_RW10434_WT_20110502_1_s3_emb1_edited",'ZD_RW10425_WT_20100412_2_s1_emb1_edited','ZD_RW10425_WT_20100412_2_s1_emb2_edited','ZD_RW10425_WT_20100412_2_s2_emb1_edited','ZD_RW10425_WT_20100412_2_s2_emb2_edited') 
wt_ashape<-dynamic_ashape(emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],is_cell=F) # a list of list of ashape 
wt_ashape_cell<-dynamic_ashape(emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],is_cell=T) # a list of list of matrices
# view in 3D
view_emb_ashape(wt_ashape[[12]],clone=1:13)
# view all time point
sapply(wt_ashape,view_emb_ashape,clone=1:13,is_save=T,dir='wt_sample1')
# To gif by ImageMagick, install by macports
cd /Volumes/baolab/xuy/distance_mx/dynamic_boundary/plot/wt_sample1
convert -delay 30 *.png -loop 0 wt_sample1.gif


# calculate connectedness based probability and the closest cell (maybe not a good way)
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


# 2015.9.30, the dynamic boundary in wwp-1 embryo 
si_info<-read.table('/Volumes/baolab/xuy/embryos_data/1368-embryos-old.csv',sep=',',as.is=T,head=T) # get ID information for RNAi embyros
si_id<-si_info[,1]
names(si_id)<-sapply(si_info[,2],function(str){strsplit(str,split='.',fixed=T)[[1]][1]})
wwp1_name<-rownames(emb_info)[grep(x=rownames(emb_info),'WWP-1')]
# For E1303
wwp1_ashape<-dynamic_ashape(emb=wwp1_name[7],end_time=emb_info[wwp1_name[7],2],is_cell=F)
wwp1_ashape_cell<-dynamic_ashape(emb=wwp1_name[7],end_time=emb_info[wwp1_name[7],2],is_cell=T)
system('mkdir plot/wwp1_e1303')
view_emb_ashape(wwp1_ashape[[125]],clone=1:13,open_new=T) # plot one to adjust view
par3d(windowRect=c(0,0,700,700)) 
mouse_mx<-par3d('userMatrix') # after adjust orietantion by mouse
mapply(view_emb_ashape,emb=wwp1_ashape,time=names(wwp1_ashape),MoreArgs=list(clone=1:13,is_save=T,dir='plot/wwp1_e1303/',open_new=F,view_mx=mouse_mx))
# cd /Volumes/baolab/xuy/distance_mx/dynamic_boundary/plot/wwp1_e1303
# convert -delay 30 *.png -loop 0 wwp1_e1303.gif
# connectedness
wwp1_connect<-lapply(wwp1_ashape_cell,cal_prob_connect)

# PIE-1, E0854
pie1_name<-rownames(emb_info)[grep(x=rownames(emb_info),'PIE-1')]
pie1_ashape<-dynamic_ashape(emb=pie1_name[2],end_time=emb_info[pie1_name[2],2],is_cell=F)
pie1_ashape_cell<-dynamic_ashape(emb=pie1_name[2],end_time=emb_info[pie1_name[2],2],is_cell=T)
system('mkdir plot/pie1_e0854')
view_emb_ashape(pie1_ashape[[length(pie1_ashape)]],clone=1:13,open_new=T) # plot one to adjust view
par3d(windowRect=c(0,0,700,700)) 
mouse_mx<-par3d('userMatrix') # after adjust orietantion by mouse
mapply(view_emb_ashape,emb=pie1_ashape,time=names(pie1_ashape),MoreArgs=list(clone=1:13,is_save=T,dir='plot/pie1_e0854/',open_new=F,view_mx=mouse_mx))
# convert -delay 30 *.png -loop 0 pie1_e0854.gif

# SKN-1, E1073
skn1_name<-rownames(emb_info)[grep(x=rownames(emb_info),'SKN-1')]
skn1_ashape<-dynamic_ashape(emb=skn1_name[1],end_time=emb_info[skn1_name[1],2],is_cell=F)
skn1_ashape_cell<-dynamic_ashape(emb=skn1_name[1],end_time=emb_info[skn1_name[1],2],is_cell=T)
system('mkdir plot/skn1_e1073')
view_emb_ashape(skn1_ashape[[length(skn1_ashape)]],clone=1:13,open_new=T) # plot one to adjust view
par3d(windowRect=c(0,0,700,700)) 
mouse_mx<-par3d('userMatrix') # after adjust orietantion by mouse
mapply(view_emb_ashape,emb=skn1_ashape,time=names(skn1_ashape),MoreArgs=list(clone=1:13,is_save=T,dir='plot/skn1_e1073/',open_new=F,view_mx=mouse_mx))
# convert -delay 30 *.png -loop 0 skn1_e1073.gif

# GLP-1, E0499
glp1_name<-rownames(emb_info)[grep(x=rownames(emb_info),'GLP-1')]
glp1_ashape<-dynamic_ashape(emb=glp1_name[4],end_time=emb_info[glp1_name[4],2],is_cell=F)
glp1_ashape_cell<-dynamic_ashape(emb=glp1_name[4],end_time=emb_info[glp1_name[4],2],is_cell=T)
system('mkdir plot/glp1_e0499')
view_emb_ashape(glp1_ashape[[length(glp1_ashape)]],clone=1:13,open_new=T) # plot one to adjust view
par3d(windowRect=c(0,0,700,700)) 
mouse_mx<-par3d('userMatrix') # after adjust orietantion by mouse
mapply(view_emb_ashape,emb=glp1_ashape,time=names(glp1_ashape),MoreArgs=list(clone=1:13,is_save=T,dir='plot/glp1_e0499/',open_new=F,view_mx=mouse_mx))
# convert -delay 30 *.png -loop 0 glp1_e0499.gif


# 2015.10.1, predict cell parameter based on membrane data in "/Volumes/baolab/xuy/embryos_data/membrane_strain/ZD_BV24_JournalV_1_s1_ed_to_twitching_2.xml". Then we can determine cell interactions.
# In that embryo, t230 is the end of gastrulation.
cell_nuc<-c(rep(1.5,8),rep(2,5)) # the ratio between cell diameter and nuclear diameter. The cell radius around division must be curated. The nucleus signal is much smaller than usual!
names(cell_nuc)<-founder
boundary_contact<-function(emb,time,boundary,all_cell=F,is_boundary=T,check_boundary=F,cell_nuc_ratio=cell_nuc){
  emb_path<-"/Volumes/baolab/xuy/embryos_data/"
  time_txt<-time
  if(time<10) time_txt<-paste('0',time_txt,sep='')
  if(time<100) time_txt<-paste('0',time_txt,sep='')
  nu_tm_raw<-read.table(paste(emb_path,emb,'/CurateDiv/t',time_txt,'-nuclei',sep=''),sep=',',as.is=T,strip.white=T) 
  nu_tm_valid<-nu_tm_raw[nu_tm_raw[,2]==1&sapply(nu_tm_raw[,10],substring,0,1)!='N',]
  nu<-nu_tm_valid[,10]
  nu_an<-sapply(nu,get_ancestor) 
  nu_tm_valid[,8]<-nu_tm_valid[,8]/trans_z # Maybe I do not need to order cells
  nu_dist<-apply(nu_tm_valid[,6:8],1,cal_euclidean,data=nu_tm_valid)
  ap<-max(nu_dist)
  nu_dist<-nu_dist/max(nu_dist) # normalize
  radius_by_closest<-sapply(rownames(nu_tm_valid),function(x){ # the rownames is CHARACTER!
    nu_dist[x,order(nu_dist[x,])[2]]/2 # First estimation, the closest distance divided by two
  })
  ratio<-cell_nuc_ratio[nu_an]
  ratio[sapply(nu_an,substring,0,2)=='AB'&sapply(nu_tm_valid[,10],nchar)<=8]<-2 # AB lineage before ABxxxxxx, set the ratio as 2
  radius_by_cell<-nu_tm_valid[,9]/2*ratio/ap # Second estimation, from nucleus with a ratio
  cell_radius<-mapply(max,radius_by_cell,radius_by_closest) # take the maximum as cell radius
  radius_plus<-sapply(cell_radius,function(x){x+cell_radius})
  contact<-matrix(nu_dist<=radius_plus,nr=length(cell_radius))
  rownames(contact)<-rownames(nu_tm_valid)->colnames(contact)
  diag(contact)<-F # set self-contact to 0
  if(all_cell){
    rownames(contact)<-nu_tm_valid[,10]->colnames(contact)
    return(contact)
  }
  if(is_boundary){
    boundary_by_lineage<-lapply(founder,function(x){
      boundary_index<-boundary[[which(founder==x)]][,5]
      lineage_contact<-contact[nu_an==x,]
      if(sum(boundary_index!=2)==0){
        inner_neighbor<-''
      }else if(sum(boundary_index!=2)==1){ 
        if(sum(lineage_contact[names(boundary_index)[boundary_index!=2],])==0) inner_neighbor<-'' else inner_neighbor<-nu_tm_valid[ lineage_contact[names(boundary_index)[boundary_index!=2],] , 10]
      }else{
        inner_neighbor<-unique(unlist(apply(lineage_contact[names(boundary_index)[boundary_index!=2],],1,function(y){list(nu_tm_valid[y,10])})))
      }
#      inner_outer_ratio<-apply(lineage_contact[names(boundary_index)[boundary_index!=2],],1,function(y){neighbor_an<-sapply(nu_tm_valid[y,10],get_ancestor);return(sum(neighbor_an!=x)/length(neighbor_an)*100)})
      if(sum(boundary_index==2)==1){ # If there is only one boundary cell, there must be only one cell in total.
        outer_neighbor<-unique(nu_tm_valid[lineage_contact,10])
      }else{ 
        outer_neighbor<-unique(unlist(apply(lineage_contact[names(boundary_index)[boundary_index==2],],1,function(y){list(nu_tm_valid[y,10])})))
      }
      if(inner_neighbor[1]!='') inner_neighbor_noself<-inner_neighbor[sapply(inner_neighbor,get_ancestor)!=x] else inner_neighbor_noself<-inner_neighbor
      outer_neighbor_noself<-outer_neighbor[sapply(outer_neighbor,get_ancestor)!=x]
      if(check_boundary&(inner_neighbor[1]!='')) return(sum(!inner_neighbor_noself %in% outer_neighbor_noself)/length(outer_neighbor_noself)*100)
      if(check_boundary&(inner_neighbor[1]=='')) return(0)
      inner_neighbor_only<-inner_neighbor_noself[! inner_neighbor_noself %in% outer_neighbor_noself]
      if(length(inner_neighbor_only)==0) inner_neighbor_only<-''
      # order neighbors
      if(length(inner_neighbor_only)>1) inner_neighbor_order<-inner_neighbor_only[order(founder_index[sapply(inner_neighbor_only,get_ancestor)],inner_neighbor_only)] else inner_neighbor_order<-inner_neighbor_only
      if(length(outer_neighbor_noself)>1) outer_neighbor_order<-outer_neighbor_noself[order(founder_index[sapply(outer_neighbor_noself,get_ancestor)],outer_neighbor_noself)] else outer_neighbor_order<-outer_neighbor_noself
      neighbors<-list(inner_neighbor_order,outer_neighbor_order)
      names(neighbors)<-c('inner','outer')
      return(neighbors)
    })
    names(boundary_by_lineage)<-founder
    return(boundary_by_lineage)
  }
}

wt_neighbor<-boundary_contact(emb=wt_sample[1],time=emb_info[wt_sample[1],2],boundary=wt_ashape_cell[[120]])

# 2015.10.2 convert neighbors into ratio and set cutoff on mixture of two lineages
neighbor_to_num<-function(neighbors,cut_inner=2,print_inner=F){
  mix<-matrix(rep(0,13*13),nr=13)
  rownames(mix)<-founder->colnames(mix)
  inner_total<-sapply(neighbors,function(x){
    mark<-rep(0,13)
    names(mark)<-founder
    if(x[[1]][1]=='') return(mark)
    inner_an<-sapply(x[[1]],get_ancestor)
    inner_sum<-tapply(inner_an,inner_an,length)
    mark[names(inner_sum)]<-inner_sum
    return(mark)
  })
  outer_total<-sapply(neighbors,function(x){
    mark<-rep(0,13)
    names(mark)<-founder
    outer_an<-sapply(x[[2]],get_ancestor)
    outer_sum<-tapply(outer_an,outer_an,length)
    mark[names(outer_sum)]<-outer_sum
    return(mark)
  })
  res<-list(inner_total,outer_total)
  names(res)<-c('inner','outer')
  if(print_inner){
    for(i in 1:13){
      for(j in 1:i){
        if(inner_total[i,j]>=cut_inner&inner_total[j,i]>=cut_inner) cat("mix:",founder[i],founder[j],"\n")
      }
    }
  }
  return(res)
}

wt_neighbor_num<-neighbor_to_num(wt_neighbor)

# compare WT boundary cell number and the aggregation of the clone in time profiling
dynamic_boundary_perc<-function(ashape_cell){
  sapply(ashape_cell,function(x){
    sapply(x,function(y){
      return(sum(y[,5]==2)/(dim(y)[1])*100)
    })
  })
}
wt_boundary_perc<-dynamic_boundary_perc(wt_ashape_cell)
# plot boundary percentage with time in 'wt_sample[1]'
image_interval<-1.25
pdf('plot/wt_sample1_bou_perc.pdf')
par(cex=2,las=1,mar=c(4,4,1,1),xpd=NA)
plot((0:((dim(wt_boundary_perc)[2])-1))*image_interval,wt_boundary_perc[1,],col=founder_col[1],ylim=c(50,100),frame=F,xlab='Time to 26-cell (min)',ylab='Boundary cell percentage (%)',lwd=4,type='l')
for(i in 2:11) lines((0:((dim(wt_boundary_perc)[2])-1))*image_interval,wt_boundary_perc[i,],col=founder_col[i],lwd=4)
dev.off()
# hard to see clearly

# GLP-1 neighbors, e0499
glp1_boundary_perc<-dynamic_boundary_perc(glp1_ashape_cell[length(glp1_ashape_cell)]) # the last time point
# Boundary percentage indicates the aggregation within lineage. The aggregation in AB lineage change a lot. The lineages with ABala fate become more clustered. ABarp loses planar feature, which is adapted to ABprp.
glp1_neighbor_num<-neighbor_to_num(boundary_contact(emb=glp1_name[4],time=emb_info[glp1_name[4],2],boundary=glp1_ashape_cell[[length(glp1_ashape_cell)]]),print_inner=T)
# The inner contact indicates the mixture between two lineages. The outer contact indicates the enviroment of the lineage.

# compare the change of outer contact
compare_out_contact<-function(si,ref,cutoff){}


# Decide to compare Boudary dynamically with time scale
dynamic_contact<-function(emb,end_time,all_boundary,is_num=F,all_cell=F){
  emb_path<-"/Volumes/baolab/xuy/embryos_data/"
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
  all_neighbor<-mapply(boundary_contact,time=start_time:end_time,boundary=all_boundary,MoreArgs=list(emb=emb,all_cell=all_cell),SIMPLIFY = F)
  if(is_num){
    neighbor_mx<-lapply(all_neighbor,neighbor_to_num)
    names(neighbor_mx)<-start_time:end_time
    return(neighbor_mx)
  }
  return(all_neighbor)
}


# 2015.10.05, adjust nucleus radius during division by CurateDiv.java
mkdir /Volumes/baolab/xuy/embryos_data/ZD_RW10348_WT_20110126_2_s2_emb1_edited/CurateDiv/
javac CurateDiv.java
java BatchCurate wt_sample1.txt
wt_neighbor<-dynamic_contact(emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],all_boundary=wt_ashape_cell) # override previous one
wt_all_neighbor<-dynamic_contact(emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],all_boundary=wt_ashape_cell,all_cell=T)

# For GLP-1 e499
system(paste('mkdir //Volumes//baolab//xuy//embryos_data//',glp1_name[4],'//CurateDiv',sep=''))
java BatchCurate glp-1.txt
glp1_neighbor<-dynamic_contact(emb=glp1_name[4],end_time=emb_info[glp1_name[4],2],all_boundary=glp1_ashape_cell)

source('neighbor_function.r')
# print well-formated neighbor
sapply(founder[1:8],print_neighbor,neighbor=glp1_neighbor,gene_name='GLP-1',boundary=glp1_ashape_cell)
sapply(founder[1:11],print_neighbor,neighbor=wt_neighbor,gene_name='WT1',boundary=wt_ashape_cell)

# all patterns in all pairs of lineage in one embryo: GLP-1
pair_founder1<-founder[as.vector(t(combn(13,2)))]
names(pair_founder1)<-paste(founder[as.vector(t(combn(13,2)))],founder[as.vector(t(combn(13,2)[2:1,]))],sep='-')
glp1_pattern<-mapply(window_neighbor,lin1=pair_founder1,lin2=founder[as.vector(t(combn(13,2)[2:1,]))],MoreArgs=list(neighbor=glp1_neighbor,emb=glp1_name[4],end_time=emb_info[glp1_name[4],2],win=1,is_pattern=T,is_tiling=F,cut_nei=0.5),SIMPLIFY=F)
glp1_PATTERN<-mapply(heatmap_neighbor_pattern,pattern=glp1_pattern,lin1=pair_founder1,lin2=founder[as.vector(t(combn(13,2)[2:1,]))],MoreArgs=list(emb='glp1',is_pattern=T,is_plot=T,plot_pattern=T),SIMPLIFY=F)

# all patterns in WT1
wt_pattern<-mapply(window_neighbor,lin1=pair_founder1,lin2=founder[as.vector(t(combn(13,2)[2:1,]))],MoreArgs=list(neighbor=wt_neighbor,emb=wt_sample[1],end_time=emb_info[wt_sample[1],2],win=1,is_pattern=T,is_tiling=F,cut_nei=0.5),SIMPLIFY=F)
wt_PATTERN<-mapply(heatmap_neighbor_pattern,pattern=wt_pattern,lin1=pair_founder1,lin2=founder[as.vector(t(combn(13,2)[2:1,]))],MoreArgs=list(emb='wt1',is_pattern=T,is_plot=T,plot_pattern=T),SIMPLIFY=F)

# pattern to matrix
glp1_mx<-pattern_to_interaction(glp1_PATTERN,end=length(glp1_pattern[[1]])) # [1] 378 378
heatmap_mx(emb='glp1',glp1_mx)
wt_mx<-pattern_to_interaction(wt_PATTERN,end=length(wt_pattern[[1]]))
heatmap_mx(emb='wt',wt_mx)

save.image('dynamic_boundary.rdata')


# 2016.3.12, clonal boundary for BV24
bv24d1_ashape<-dynamic_ashape(emb='membrane_strain',end_time=230,is_cell=F)
view_emb_ashape(bv24d1_ashape[[230-73]],clone=1:13)
rgl.postscript(filename='plot/BV24_D1/t230.pdf',fmt='pdf')
