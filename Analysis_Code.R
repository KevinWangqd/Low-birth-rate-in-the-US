#library(usmap)
#library(patchwork)
library(ggplot2)
library(cleaner)
library(raster)
library(sf)
library(rgdal)
library(dplyr)
library(INLA)
library(viridis)

#Read data: d for whole dataset and excess for extremely low dataset (by u=35, trate>=35)
d <- read.csv("data_scale.csv")

#boxplot(Birth.Rate ~ Year, data = d, range = 0, xlab = "Year", ylab = "Birth rate", col = "darkseagreen", border = "navy", main = "")

#d2009 <- d[which(d$Year==2009),]
#d2019 <- d[which(d$Year==2019),]
#summary(d2009)

#summary_data <- d2019 %>%
#  summarise_all(list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE)))

# Print the summary dataframe
#print(summary_data)



#a <- as.data.frame(unique(excess$index))
#a$location <- sub(", [0-9]+$", "", a$`unique(excess$index)`)
#b <- freq(a$location)
#b1 <- b[b$count>10,]
#write.csv(b1,"index_10year.csv")


u <- 35
#u <- 38




excess <- d[which(d$trate>=u),]
#countyyear <- read.csv("index_10year.csv")
#excess <- excess[excess$County %in% countyyear$item,]
#excess <- read.csv("excess_scale.csv")
#u=13
excess$rate <- excess$trate-u

excess$Crisis <- ifelse(excess$Year<=2008,0,1)
excess$Crisis <- as.factor(excess$Crisis)





#excess$Mother.s.Hispanic.Origin <- as.factor(excess$Mother.s.Hispanic.Origin)

#get US map
m <- getData(name = "GADM", country = "USA", level = 0)

m <- m %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(.)) %>%
  arrange(desc(area)) %>%
  slice(1)

plot(m)
## Mesh construction
coo <- cbind(excess$long,excess$lat)
max.edge <- diff(range(excess$long))/15 ##max.edge value to be between 1/3 to 1/10 times smaller than the spatial range
bound.outer <- diff(range(excess$long))/3 ## 1/3 of spatial range
mesh <- inla.mesh.2d(loc = coo,
                     max.edge = c(1,2)*max.edge, # -use 5 times max.edge in the outer extension/offset/boundary
                     cutoff = max.edge/5,
                     offset = c(max.edge,bound.outer))
plot(mesh)
points(coo,col="red")
mesh$n

## Index: train and validation
index <- excess$Year
train_index <- which(excess$Year <= 2019)
test_index <- which(excess$Year==2020)

## SPDE (PC prior)
spde <- inla.spde2.pcmatern(mesh=mesh,
                            alpha = 2,
                            prior.range = c(10,0.01), #P(range<10=0.01
                            prior.sigma = c(1,0.01) # P(sigma>1)=0.01
)


timesn <- length(unique(excess$Year))

indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn
)
lengths(indexs)

## Projection matrix A (training) and Ap (prediction/validation)

group <- excess$Year - min(excess$Year) + 1
A <- inla.spde.make.A(mesh = mesh, loc = coo[train_index, ], group = group[train_index], n.group = timesn)
Ap <- inla.spde.make.A(mesh = mesh, loc = coo[test_index, ], group = group[test_index], n.group = timesn)
dim(A)


## Stack (estimation and prediction/validation)
stk.e <- inla.stack(
  tag = "train",
  data = list(y=excess$rate[train_index]),
  A = list(1, A),
  effects = list(data.frame(Intercept = rep(1, length(train_index)),
                            Year=excess$Year[train_index],
                            race=excess$Mother.s.Hispanic.Origin[train_index],
                            id=excess$X[train_index],
                            Long=excess$long_s[train_index],
                            Lat=excess$lat_s[train_index],
                            birth_weight=excess$Average.Birth.Weight[train_index],
                            Precipitation=excess$Precipitation[train_index],
                            Temperature=excess$Temperature[train_index],
                            GDP=excess$GDP_per[train_index],
                            PI=excess$PI[train_index],
                            RUCC=excess$RUCC[train_index],
                            Marriage=excess$Crude.Marriage.Rate[train_index],
                            Elder=excess$Elder_pop[train_index],
                            Edu=excess$bachelor_percentage[train_index],
                            Female_proportion=excess$female.proportion[train_index],
                            Age_of_Mother=excess$Average.Age.of.Mother[train_index],
                            LMP_Gestational_Age=excess$Average.LMP.Gestational.Age[train_index]), s = indexs)
)


stk.p <- inla.stack(
  tag = "test",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(Intercept = rep(1, length(test_index)),
                            Year=excess$Year[test_index],
                            race=excess$Mother.s.Hispanic.Origin[test_index],
                            id=excess$X[test_index],
                            Long=excess$long_s[test_index],
                            Lat=excess$lat_s[test_index],
                            birth_weight=excess$Average.Birth.Weight[test_index],
                            Precipitation=excess$Precipitation[test_index],
                            Temperature=excess$Temperature[test_index],
                            GDP=excess$GDP_per[test_index],
                            PI=excess$PI[test_index],
                            RUCC=excess$RUCC[test_index],
                            Marriage=excess$Crude.Marriage.Rate[test_index],
                            Elder=excess$Elder_pop[test_index],
                            Edu=excess$bachelor_percentage[test_index],
                            Female_proportion=excess$female.proportion[test_index],
                            Age_of_Mother=excess$Average.Age.of.Mother[test_index],
                            LMP_Gestational_Age=excess$Average.LMP.Gestational.Age[test_index]), s = indexs)
)

stk.full <- inla.stack(stk.e, stk.p)


## Prior
hyper = list(theta = list(prior="pc.gevtail", param=c(7, 0, 0.8)))

rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))

## formula and run INLA
fromula1 <- y ~ -1+Intercept + race + Lat +Long + birth_weight  + RUCC + Elder + Edu +
  Precipitation + Temperature + GDP + PI + Female_proportion + Age_of_Mother + Marriage+
  LMP_Gestational_Age + f(id, model = "iid") 

#fromula2 <- y ~ -1 + Intercept + race + Lat +Long + birth_weight  + RUCC + Elder + Edu +
#  Precipitation + Temperature + GDP + PI + Female_proportion + Age_of_Mother + Marriage+
#  LMP_Gestational_Age + f(Year, model="ar1") 


fromula4 <- y ~ -1+Intercept + race + Lat +Long + birth_weight  + RUCC + Elder + Edu +
  Precipitation + Temperature + GDP + PI + Female_proportion + Age_of_Mother + Marriage+
  LMP_Gestational_Age + f(s, model = spde, group = s.group, control.group = list(model = "ar1", hyper = rprior))+ f(id, model = "iid")

fromula5 <- y ~ -1+Intercept + race + Lat +Long + birth_weight  + RUCC + Elder + Edu +
  Precipitation + Temperature + GDP + PI + Female_proportion + Age_of_Mother + Marriage+
  LMP_Gestational_Age + f(Year, model="ar1") + f(s, model = spde)+ f(id, model = "iid") 


result <- inla(fromula4, 
               data = inla.stack.data(stk.full),
               family= "gp", 
               control.family = list(list(control.link  = list(model = "quantile", quantile = 0.65)#, hyper = list(theta  = hyper.pre)
                                          )), 
               control.predictor=list(A=inla.stack.A(stk.full), link=1, compute=TRUE),  
               control.compute = list(dic = TRUE, waic = TRUE, config = TRUE, cpo = TRUE),
              verbose = T)



summary(result)


result$summary.hyperpar

hist(result$cpo$pit)


## Result processing and data preparation for evaluation
index_train_inla <- inla.stack.index(stk.full, "train")$data
index_val_inla <- inla.stack.index(stk.full, "test")$data

result.train1 <- result$summary.fitted$mean[index_train_inla]
result.val1 <- result$summary.fitted$mean[index_val_inla]

M_fit_spde1 <- array(NA, nrow(excess))
M_fit_spde1[train_index] <- result.train1
M_fit_spde1[test_index] <- result.val1


# Model Comparison (We only include the code for Model 1. For Model 2/3/4, the user may change the likelihood and the formula)

## DIC and WAIC
result$dic$dic
result$waic$waic

## RMSE
## Functions (slcpo is the same as Logarithm score)
RMSE <- function(set, outcome, data, fit) {
  res <- data[set, outcome] - fit[set]
  RMSE_val <- sqrt(mean(res^2, na.rm = TRUE))
  return(RMSE_val)
}


rmse_train_spde1 <- RMSE(train_index, "rate", excess, M_fit_spde1)
rmse_val_spde1 <- RMSE(test_index, "rate", excess, M_fit_spde1)
rmse_train_spde1
rmse_val_spde1

## CPO and PIT
hist(result$cpo$pit, xlab = "", main = "")

slcpo <- function(m) {
  -sum(log(m$cpo$cpo),na.rm=TRUE)
}

cpo_1 <- slcpo(result)
cpo_1
##visualization

plotting_data1 <- data.frame(Observed = excess$rate,
                             Predicted = M_fit_spde1)
plotting_data1$group[train_index] <- "train"
plotting_data1$group[test_index] <- "validation" 
plotting_data1$group[train_index] <- "train"





p1 <- ggplot(data = plotting_data1, aes(x = Observed, y = Predicted)) +
  geom_point(aes(color = group)) + geom_abline(intercept = 0, slope = 1) + xlim(0, 10) +
  ylim(0, 10) + theme_bw() + ggtitle("")
p1


## Correlation
obs_val <- plotting_data1[which(plotting_data1$group=="validation"),1]
pre_val <- plotting_data1[which(plotting_data1$group=="validation"),2]
scatter.smooth(obs_val,pre_val)
cor(obs_val,pre_val) 



## The following is about the plot for spatial effect

bb <- st_bbox(m)
x <- seq(bb$xmin - 1, bb$xmax + 1, length.out = 200)
y <- seq(bb$ymin - 1, bb$ymax + 1, length.out = 200)
dp <- as.matrix(expand.grid(x, y))
plot(dp, asp = 1)
p <- st_as_sf(data.frame(x = dp[, 1], y = dp[, 2]),
              coords = c("x", "y")
)
st_crs(p) <- st_crs(4326)
ind <- st_intersects(m,p)

dp <- dp[ind[[1]], ]
plot(dp, asp = 1)

dp <- as.data.frame(dp)

write.csv(dp,"dp.csv")


spatial <- read.csv("dp.csv")
coord <- spatial[,2:3]
coord <- as.matrix(coord)
coord <- SpatialPoints(coord, proj4string = CRS(as.character(NA)),
                       bbox = NULL)

gproj <- inla.mesh.projector(mesh,  coord)

g.mean <- inla.mesh.project(gproj, result$summary.random$s$mean[1:mesh$n], projection="longlat")
g.sd <- inla.mesh.project(gproj, result$summary.random$s$sd[1:mesh$n])


a <- as.data.frame(g.mean)
b <- as.data.frame(g.sd)

spatial$gmean <- a$g.mean
spatial$gsd <- b$g.sd

p1 <- ggplot(m) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = spatial, aes(x = long, y = lat, color = gmean),
    size = 2
  )  +
  labs(x = "", y = "",colour = "") +
  scale_color_viridis(option="inferno",begin = 0.2) +  theme_bw()+ggtitle("")





p2 <- ggplot(m) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = spatial, aes(x = long, y = lat, color = gsd),
    size = 2
  ) + 
  labs(x = "", y = "",colour = "") +scale_color_viridis(option="inferno",begin = 0.2)+theme_bw()+ggtitle("")





#jpeg("ModelV.jpeg", quality = 100, units = "in", width = 6, height = 6, res = 300)

# Print the plot or visualization
#p1
# Close the jpeg device
#dev.off()






##Stwscrps


library("evd")

## Stwcrps for GP dist. note here the y is excess

stwcrps_gp = function(y, q, ξ, p) {
  S = abs(expected_twcrps_gp(q, ξ, p))
  twcrps = twcrps_gp(y, q,ξ, p)
  twcrps / S + log(S)
}


twcrps_gp = function(y, q,  ξ, p, q_level=0.65) {
  ## This is only for the reparametrization of σ and q (50% quantile)
  #q_level=q_level
  if (ξ==0) {
    σ = -q/log(1-q_level)
  } else {
    σ = q*ξ/((1-q_level)^(-ξ)-1)}
  
  F = function(x) sapply(x, function(z) mean(pgpd(z, loc=0, scale=σ, shape=ξ)))
  quantiles = sapply(σ, function(σ) qgpd(p, loc=0, scale=σ, shape=ξ))
  
  if (length(quantiles) == 1) {
    y_min = quantiles
  } else {
    y_min = uniroot(function(x) F(x) - p, lower = min(quantiles), upper = max(quantiles))$root
  }
  p_max = .999
  y_max = max(qgpd(p_max, loc=0, scale=σ, shape=ξ), max(y) + 1)
  res = rep(0, length(y))
  for (i in seq_along(y)) {
    if (y[i] < y_min) {
      res[i] = integrate(function(x) (1 - F(x))^2, lower = y_min, upper = y_max)$value
    } else if (y[i] < y_max) {
      res[i] = integrate(function(x) F(x)^2, lower = y_min, upper = y[i])$value
      res[i] = res[i] + integrate(function(x) (1 - F(x))^2, lower = y[i], upper = y_max)$value
    } else {
      res[i] = integrate(function(x) F(x)^2, lower = y_min, upper = y_max)$value
    }
  }
  res = res + (y_min - y) * (ifelse(y <= y_min, 1, 0) - p)^2
  res
}





expected_twcrps_gp = function(q, ξ, p, 
                              q_true = q, ξ_true = ξ, q_level=0.65) {
  σ_true = ifelse(ξ_true==0,-q_true/log(1-q_level), q_true*ξ_true/((1-q_level)^(-ξ_true)-1))
  
  p_min = .00001
  y_min = min(qgpd(p_min,loc=0, scale=σ_true, shape=ξ_true))
  p_max = .99999
  y_max = max(qgpd(p_max,loc=0, scale=σ_true, shape=ξ_true))
  
  if (length(c(q_true, ξ_true)) == 2) {
    density = function(x) dgpd(x, loc=0, scale=σ_true, shape=ξ_true)
  } else {
    density = function(x) sapply(x, function(z) mean(dgpd(z, loc=0, scale= σ_true, shape= ξ_true)))
  }
  integrate(function(y) density(y) * twcrps_gp(y, q, ξ, p),
            lower = y_min, upper = y_max)$value
}



stwcrps_train= mean(stwcrps_gp(y=excess$rate[index_train_inla],q=result.train1,ξ=result$summary.hyperpar$mean[1], p=0.7))
stwcrps_train

stwcrps_val=mean(stwcrps_gp(y=excess$rate[index_val_inla],q=result.val1,ξ=result$summary.hyperpar$mean[1],p=0.7))
stwcrps_val


#stwcrps_train=mean(stwcrps_gp(y=excess$rate[index_train_inla],q=result.train1,ξ=result$summary.hyperpar$mean[1],p=0.8))
#stwcrps_train

#stwcrps_val=mean(stwcrps_gp(y=excess$rate[index_val_inla],q=result.val1,ξ=result$summary.hyperpar$mean[1],p=0.8))
#stwcrps_val






















## Appedix POT select
library(POT)



rate <- na.omit(d$trate)

mrplot_1 <- mrlplot(rate)
data <- data.frame(mrplot_1)
max(mrplot_1$x)

mrlplot(rate,u.range = c(30,42.56))

par(mfrow = c(1,1))

tc <- POT::tcplot(rate, nt = 150)

plot(tc[[1]])
plot(tc[[2]])


tc <- POT::tcplot(d$trate,nt=50,)
# Assuming you already have the tcplot stored in the variable 'tc'

plot(tc[[1]])

# Plot the second panel
plot(tc[[2]])
# Add vertical lines at x = 3.67
abline(v = 3.67, col = "red", lty = 2)  # Add a red dashed vertical line

# Create the mrlplot
mrplot_1 <- mrlplot(National$T90)

# Add the custom line with slope -0.3 and intercept 5.5
mrlplot(rate,u.range = c(30,42.56))
abline(a = 18.65, b = -0.47, col = "blue", lty = 2)
abline(v = 35, col = "red", lty = 2)  # Add a red dashed vertical line


jpeg("mrl1.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)

mrlplot(rate,u.range = c(25,42.56))
abline(a = 18.65, b = -0.47, col = "blue", lty = 2)
abline(v = 35, col = "red", lty = 2)  # Add a red dashed vertical line

dev.off()



jpeg("mrl_scale.jpeg", quality = 100, units = "in", width = 6, height = 6, res = 300)

# Print the plot or visualization
POT::tcplot(rate, nt=100,lty = 1, u.range = c(30,43),which = 1)

# Add the custom line with slope -0.3 and intercept 5.5
abline(v = 35, col = "red", lty = 2)  # Add a red dashed vertical line
# Close the jpeg device
dev.off()







