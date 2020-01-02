#All values in the original data should be real numbers according to real-life experiences.

#Solution to problem 3a.
#Input:traj--a n*3 matrix,the first two columns stand for the position while the third stands for time. 
#Output:traj_trans--the former matrix with time starting from 0 and position starting at the origin.
translate=function(traj){
  x_start=traj[1,1]
  y_start=traj[1,2]
  t_start=traj[1,3]
  traj_trans=traj
  starts=c(x_start,y_start,t_start)
  for(i in 1:nrow(traj)){
    for(j in 1:3){
      traj_trans[i,j]=traj[i,j]-starts[j]
    }
  }
  return(traj_trans)
}

#Solution to problem 3b.
#Input:traj--a n*3 matrix with the first row being (0,0,0).
#Output:theta--the angle (in rads) formed by the secant line connecting the origin and the final position;values between -pi and pi.
angle_compute=function(traj){
  final=nrow(traj)
  theta=acos(traj[final,1]/sqrt(traj[final,1]^2+traj[final,2]^2))*sign(traj[final,2])
  return(theta)
}

#Solution to problem 3c.
#Input:traj--a n*3 matrix with the first row being (0,0,0).
#Output:traj_rotate--the rotated matrix with its final point on the positive x-axis.
#Definition from problem 3b is used.
#The y-value of the final spot might be a minute number due to approximate calculations.
#So you may use round(rotate(input),n) to check the results, where n is a (rather small) positive integer. 
rotate=function(traj){
  theta=angle_compute(traj)
  traj_rotate=traj
  for(i in 1:nrow(traj)){
    traj_rotate[i,1]=traj[i,1]*cos(theta)+traj[i,2]*sin(theta)
    traj_rotate[i,2]=traj[i,2]*cos(theta)-traj[i,1]*sin(theta)
  }
  return(traj_rotate)
}

#Solution to problem 3d.
#Input:traj--a n*3 matrix.
#Output:traj_norm--the normalized matrix which begins at the origin at time 0 and ends on the positive x-axis.
#Definitions from the previous problems are used.
normalize=function(traj){
  traj_trans=translate(traj)
  traj_norm=rotate(traj_trans)  
  return(traj_norm)
}

#Solutions to problem 3e.
#Input:traj--a normalized n*3 matrix.
#Output:dist--the total (Euclidean) distance traveled.
#       max_dev--the maximum absolute deviation from the secant.
#       ave_dev--the average absolute deviation of the observed trajectory from the direct path.
#       area--the absolute area under the curve for the trajectory relative to the secant line.
metric_calc=function(traj){
  index=nrow(traj)
  dist_step=rep(0,index-1)
  for(i in 2:index){
    dist_step[i-1]=sqrt((traj[i,1]-traj[i-1,1])^2+(traj[i,2]-traj[i-1,2])^2)
  }
  dist=sum(dist_step)
  max_dev=max(abs(traj[,2]))
  ave_dev=mean(abs(traj[,2]))
  area_step=rep(0,index-1)
  for(i in 2:index){
    area_step[i-1]=(abs(traj[i,2])+abs(traj[i-1,2]))*(traj[i,1]-traj[i-1,1])/2
  }
  area=sum(area_step)
  result_list=c(dist,max_dev,ave_dev,area)
  return(result_list) 
}

#Getting results for the train set and compare with true values.
#Since the number of digits of the true values are uncertain (and unequal),
#both results are rounded to the tenth digit before comparing.
#This may introduce some error but is GENERALLY helpful in checking.
train=read.csv("train_trajectories.csv",header=TRUE)
tres=read.csv("train_measures.csv",header=TRUE)
subs=max(train[,1])
trials=max(train[,2])
eq=matrix(NA,nrow=trials,ncol=subs)
for(i in 1:subs){
  for(j in 1:trials){
    tset=train[train[,1]==i&train[,2]==j,3:5]
    if(dim(tset)[1]==0) next
    print(paste("Testing subject",i,",trial",j,"."))
    train_res=metric_calc(normalize(tset));print(train_res)
    eq[j,i]=all(round(tres[tres[,1]==i&tres[,2]==j,3:6],1)==round(train_res,1))
  }
}

#Solution to problem 3f.
#Input:N/A
#Output:metrics--a n*4 matrix concerning the traits in problem 3e, 
#where n stands for the total number of trials.
#"Intersect" is used to get a set of the subjects and trials,
#as I haven't found a better command for that.
test=read.csv("test_trajectories.csv",header=TRUE)
test_subs=intersect(test[,1],test[,1])
test_trials=intersect(test[,2],test[,2])
a=length(test_subs);b=length(test_trials);n=a*b
tlist=matrix(NA,nrow=a,ncol=b)
for(i in 1:a){
  for(j in 1:b){
    tlist[i,j]=paste("Sub",test_subs[i],"Trial",test_trials[j])
  }
}
metrics=matrix(NA,nrow=length(tlist),ncol=4,dimnames=list(tlist,c("dist","max_dev","ave_dev","area")))
for(i in 1:a){
  for(j in 1:b){
    testset=test[test[,1]==test_subs[i]&test[,2]==test_trials[j],3:5]
    if(dim(testset)[1]==0) next
    print(paste("Testing subject",test_subs[i],",trial",test_trials[j],"."))
    test_res=metric_calc(normalize(testset));print(test_res)
    metrics[b*(i-1)+j,]=test_res
  }
}
(metrics=metrics[complete.cases(metrics),])

