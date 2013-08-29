#This code allows users to take the output from morphological clines fitted in CFit 7 (Gay et al 2008*) and plot allele frequency and morphological clines, including the ability to make heat map style representations of the trait variance expected under the fitted cline.

#*Gay et al 2008 Comparing clines on molecular and phenotypic traits in hybrid zones: a window on tension zone models. Evolution 62:2789)


# Input data (for example plots)

#the data are four morphological traits through a hybrid zone between stream and anadromous sticklebacks in Bonsall Creek, Vancouver Island. The anadromous population is to the left in the graph, the x axis is centred 2.3 km from the sea 

cfitdata<-read.table(file.choose(), header=T)  #choose cfitdata.txt
attach(cfitdata)


#****plotting unimodal clines****

#p plots just the allele frequency cline as a function of x

p<-function(x,cen,w) exp(w*(x-cen))/(1+exp(w*(x-cen)))

#here's the cline data for Eda (a gene controlling stickleback armour) in the hybrid zone

centreEda<-0.008;widthEda<--2.98;

#here's plot of this cline (I just plot the points in white and then use points(type="l") to make it a smooth curve); I'm sure you can do this more elegantly.

dist<-seq(-0.8,2,by=0.01)

plot(dist,p(dist,centreEda,widthEda), col="white", axes=T, xlab="distance from sea (km)", ylab="freq of marine allele/trait")
lines(dist,p(dist,centreEda,widthEda),lty=1,col="blue")


#moy plots the a morphological trait average cline as a function of x

moy<-function(x,cen,w,min,delta) delta*p(x,cen,w)+min

#varp describes the expected variance around the cline as a function of x. s1, s2 and s3 are the variances for the left, centre and right hand parts of the cline, respectively.

varp<-function(x,cen,w,s1,s2,s3) (s1*p(x,cen,w)^2)+(2*s2*p(x,cen,w)*(1-p(x,cen,w)))+(s3*(1-p(x,cen,w))^2)

#unimodal generates a probability density for the fitted cline (the lighter the color, the more likely you are to observe datapoints at that point). min and delta are the trait minimum and the amount it changes across the cline, respectively.

unimodal<-function(x,y,cen,w,min,delta,s1,s2,s3) (exp(-((y-moy(x,cen,w,min,delta))^2)/(2*(varp(x,cen,w,s1,s2,s3)^2))))/(sqrt(2*pi)*varp(x,cen,w,s1,s2,s3))


#in the unimodal.array function, xmin, xmax, ymin ymax set the range of values for which you wish to generate the plot, and the remaining input comes from CFit 7. i.e.

# cen w 0 0 0 0
# min delta s1 s2 s3 0 0
# 0 0 0 0 0


unimodal.array<-function(xmin,xmax,ymin,ymax,cen,w,min,delta,s1,s2,s3) {
	xx<-seq(xmin,xmax,(xmax-xmin)/50)
	yy<-seq(ymin,ymax,(ymax-ymin)/50)
	
	arr<-array(0,c(length(xx),length(yy)))
	
	for(i in 1:length(xx)) {
		for(j in 1:length(yy)) 
				{arr[i,j]<-unimodal(xx[i],yy[j],cen,w,min,delta,s1,s2,s3)
				}}
				arr
				}

#this is the best fit data for the keel width cline

keelwidthcline<-unimodal.array(-1,2,-1.6,1.6,0.32,-33.56,-0.19,0.31,0.374,0.374,0.374)

#this plots the contours and the data 
quartz()
filled.contour(
		x = seq(-1,2,length.out=nrow(keelwidthcline)),
		y = seq(-1.6,1.6,length.out=ncol(keelwidthcline)),
		keelwidthcline, color = heat.colors, plot.title = title(main = "keel width cline"),nlevels=10,plot.axes={ axis(1); axis(2); points(site/1000,keelres)})




#****trimodal plotting functions****

#these functions allow for a three part cline (the below are called by trimodal.array)

pentefreqfunc<-function(slopetrait,slopemoymix,slopemoy) (4*slopetrait*slopemoy*exp(slopetrait*slopemoymix)-((1+exp(slopetrait*slopemoymix))^2)*slopetrait)/(-1+exp(2*slopetrait*slopemoymix))

f<-function(x,centre,pentefreq,xpos1,xpos2,tslope1,tslope2) {		if(x<=(centre-xpos2)) {1-((1-exp(-pentefreq*xpos2))/(1+exp(-pentefreq*xpos2)))*exp(-(tslope2*pentefreq)/(1+exp(xpos2*pentefreq))*(x-centre+xpos2))} 
		if(x>=(centre+xpos1)) {exp(pentefreq*xpos1)/(1+exp(pentefreq*xpos1))*exp(-(-(tslope1*pentefreq)/(1+exp(xpos1*pentefreq)))*(x-centre-xpos1))}
			else {exp(pentefreq*(x-centre))/(1+exp(pentefreq*(x-centre)))}
			}
			
m1<-function(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy) {f(x-slopemoymix,centre,slopemoy*slopetrait,1000,1000,tslope1,tslope2)*delta+min}

m3<-function(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy) {f(x+slopemoymix,centre,slopemoy*slopetrait,1000,1000,tslope1,tslope2)*delta+min}

m2<-function(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy) {(m1(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy)+m3(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy))/2}

ff<-function(p,f1,f2,f3,f4,f5) {
	if(p<(1/6)) {6*p*f1}
	if(p<(1/3)) {6*(f1-f2)*p+2*f1-f2}
	if(p<(1/2)) {6*(f3-f2)*p+3*f2-2*f3}
	if(p<(2/3)) {6*(f4-f3)*p+4*f3-3*f4}
	if(p<(5/6)) {6*(f5-f4)*p+5*f4-4*f5}
	else {6*(1-p)*f5}
	}

trimodal<-function(x,y,centre,slopetrait,xpos1,xpos2,tslope1,tslope2,min,delta,s1,s2,s3,slopemoymix,slopemoy,f1,f2,f3,f4,f5) {
	p<-f(x,centre,pentefreqfunc(slopetrait,slopemoymix,slopemoy),xpos1,xpos2,tslope1,tslope2)
	p1<-p^2+p*(1-p)*ff(p,f1,f2,f3,f4,f5)
	p3<-(1-p)^2+p*(1-p)*ff(p,f1,f2,f3,f4,f5)
	
	p1*(exp(-(y-m1(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy))^2/((2*s1)^2))/(sqrt(2*pi*s1)))+p3*(exp(-(y-m3(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy))^2/((2*s3)^2))/(sqrt(2*pi*s3)))+(1-p1-p3)*(exp(-(y-m2(x,centre,slopetrait,tslope1,tslope2,min,delta,slopemoymix,slopemoy))^2/((2*s2)^2))/(sqrt(2*pi*s2)))
	}


#this is the only function you need to run

trimodal.array<-function(xmin,xmax,ymin,ymax,var_list) {
	xx<-seq(xmin,xmax,(xmax-xmin)/50)
	yy<-seq(ymin,ymax,(ymax-ymin)/50)
	
	arr<-array(0,c(length(xx),length(yy)))
	
	for(i in 1:length(xx)) {
		for(j in 1:length(yy)) 
				{arr[i,j]<-trimodal(xx[i],yy[j],var_list[1],var_list[2],var_list[3],var_list[4],var_list[5],var_list[6],var_list[7],var_list[8],var_list[9],var_list[10],var_list[11],var_list[12],var_list[13],var_list[14],var_list[15],var_list[16],var_list[17],var_list[18])
				}}
				arr
				}
				
				
				
				
#trial finres data (the length of the dorsal fin through the HZ)
#the data come from the CFit maximum model:
# centre width xpos1 xpos2 tslope1 tslope2 coeffslopemean
# min delta sigma1 sigma2 sigma3 mix
# phi1 phi2 phi3 phi4 phi5

#I recommend entering your values as a list like this (the example below is for my trait finres) and using trimodal.array() to plot the clines. The variables are in the order a-r as in the table on page 10 of the CFit 7 documents,it's important to maintain this order.


finresdata<-c(centre<--0.34,slope<--2.04,d1<-1.17,d2<-0.0001,t1<-0.250,t2<-0.0001,min_mu<--0.630,delta_mu<-2.254,sigma1<-1.596,sigma2<-0.197,sigma3<-0.909,mix<-0.726,a_m<-0.555,f1<-0.968,f2<-0.882,f3<--0.581,f4<--0.013,f5<-0.675)

#the finres full model gives a really strange result

finresplot<-trimodal.array(-0.8,2,-3,3,finresdata)

xlimits<-seq(-0.8,2,length.out=nrow(finresplot))
ylimits<-seq(-3,3,length.out=ncol(finresplot))

filled.contour(xlimits,ylimits,finresplot, color = heat.colors, plot.title = title(main = "dorsal fin length"),nlevels=10,plot.axes={ axis(1); axis(2); points(site/1000,finres)})


#trial pectres data (the length of the pectoral fin through the HZ)

pectresdata<-c(centre_pectres<--0.172321142,pente_pectres<--35.19,xpos1_pectres<-0.0001,xpos2_pectres<-43.506,tslope1_pectres<-0.000,tslope2_pectres<-0.765,coeffslopemoy_pectres<-0.291,min_pectres<--0.742,delta_pectres<-1.327,sigma1_pectres<-0.737,sigma2_pectres<-0.545,sigma3_pectres<-1.698,mixpectres<-1.284,f1_pectres<-0.817,f2_pectres<--0.499,f3_pectres<--0.864,f4_pectres<--0.497,f5_pectres<--0.188)


#the pectres full model gives a two step cline (which happens because it tries to accommodate the strangely high points at 1.7)

pectresplot<-trimodal.array(-0.8,2,-3,3,pectresdata)

xlimits<-seq(-0.8,2,length.out=nrow(pectresplot))
ylimits<-seq(-3,3,length.out=ncol(pectresplot))

test.plot<-filled.contour(xlimits,ylimits,pectresplot, color = heat.colors, plot.title = title(main = "pectoral fin length"),nlevels=10,plot.axes={ axis(1); axis(2); points(site/1000,pectres)})
