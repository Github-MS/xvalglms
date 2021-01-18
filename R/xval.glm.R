#' Run X-val on a set of glms
#'
#' @param data data.frame containing the data
#' @param models list of model formulas, each list entry is used as a
#' @param glm.family name of the glm family (see \link{glm}() for details), default is gaussian
#' @param folds number of folds (default = 10)
#' @param repeats number of repeats (default = 200)
#' @param loss loss function for the GLM (default = NULL equals RMSE)
#' @param numCore number of cores for use with parallel (default = NULL, equals no parallelization)
#' @param plots output fancy plots with results
#' @param gray output greyscale plots (default = F)
#' @param seed seed for the folds (default = NULL, seed will be random)
#' @param showConsoleOutput show console output, set to FALSE to suppress all conole output (default = TRUE)
#'
#'
#' @return object of class \code{xval.glm}
#'
#' @export
xval.glm <- function(data, models, glm.family = gaussian, folds = 10, repeats = 200, loss = NULL, numCore = NULL, plots = T, gray = F, seed = NULL, showConsoleOutput = T) {
  # -------------------------------------------------------------
  # function to do K-fold cross validation on a set of glms
  # data is a data frame
  # models is a list with several model formulae
  # family is the model family used in glm()
  # K determines the number of folds
  # loss determines the LOSS function (set to standard if NULL)
  # repeats gives the number of repetitions of cross validation loops
  # @ M. de Rooij, 17-11-2017 (original code)
  # @ W.D. Weeda, 09-04-2018 (parallelization code)
  # --------------------------------------------------------------

  #necessary checks (data, models and glm.family)
  if(!is.data.frame(data)) stop('data must be a data.frame')

  if(is.list(models)) {
    if(!all(unlist(lapply(list(a~1,a~2),class))=='formula')) {
      stop('Supplied models list must contain formula')
    }
  } else {
    stop('`models` must be a list of R formula')
  }

  if(class(glm.family())!='family') stop('`glm.family` error, please provide a valid link function for glm()')

  #perform sanity checks on data (give pointers on folds)
  if(folds > nrow(data)) stop(paste0('Number of folds (',folds,') is larger than available datapoints (',nrow(data),')'))

  #set number of cores
  if(!is.null(numCore))  {

    numCore <- as.numeric(numCore)
    maxCore <- detectCores()

    if(numCore > maxCore) stop('Number of cores is higher than available.')

    cl <- makeCluster(numCore)

    if(showConsoleOutput) cat('Using',numCore,'cores.\n')

    #register cluster
    registerDoParallel(cl)
  }

  #set seeds
  if(is.null(seed)) seed <- runif(1,0,1000)
  set.seed(seed)

  #define loss label for plot
  if(is.null(loss)) my.ylab <- "RMSE_P" else  my.ylab = "PE"

  #define lengths of loops
  K <- folds
  M <- length(models)
  n <- nrow(data)
  y <- as.matrix(model.response(model.frame(models[[1]], data = data)))
  folds <- cut(seq(1,n),breaks = K,labels = FALSE)

  #start time
  tstart <- Sys.time()

  #make output lists
  out <- list()
  preds <- array(NA,dim=c(n,repeats,M))

  #parallel loop of repeats
  if(!is.null(numCore)) {

    if(showConsoleOutput) cat('Running Cross-validation...')

    tot_cval_out <- foreach(it=1:repeats,.combine = c) %dopar% {
    set.seed(it+seed)
    cval_out <- cross.validate(M, folds, n, K, glm.family, data, y, models, loss)

      return(list(cval_out))
    }

    #get parallel output in correct output objects
    for(i in 1:repeats) {
      out <- c(out,tot_cval_out[[i]]$loss)
      preds[,i,] <- tot_cval_out[[i]]$pred
    }

    #stop parallel loop
    stopCluster(cl)

  } else {

    if(showConsoleOutput) cat('Running Cross-validation...\n')
    if(showConsoleOutput) pbar <- txtProgressBar(1,repeats,1,style=3)

    for(i in 1:repeats) {
      cval_out <- cross.validate(M, folds, n, K, glm.family, data, y, models, loss)
      out <- c(out,cval_out$loss)
      preds[,i,] <- cval_out$pred
      if(showConsoleOutput) setTxtProgressBar(pbar,i)
    }
    if(showConsoleOutput) cat('\n')
  }

  #stop time
  tend <- Sys.time()
  if(showConsoleOutput) cat(paste0('done [ ',round(as.numeric(difftime(tend,tstart,units='sec')),1),' sec ]\n'))

  #put all repeats in nice matrix
  RMSEP <- matrix(unlist(out), ncol = M, nrow = repeats, byrow = T)

  #count the number of wins of each model over the repeats
  winners <- apply(RMSEP,1,which.min)
  winmat <- matrix(0, ncol = M, nrow = repeats)
  for(i in 1:length(winners)){
    winmat[i,winners[i]] <- 1
  }
  wins <- apply(winmat,2,sum)

  #concatenate cross-validation results
  cv.pe <- matrix(NA,M*repeats,2)
  colnames(cv.pe) <- c("Model", "RMSEP")
  cv.pe[,1] <- rep(seq(1,M), each = repeats)
  cv.pe[,2] <- RMSEP
  cv.pe <- as.data.frame(cv.pe)

  #check stability (cumulative proportion of wins over repeats)
  stab <- data.frame()
  startstab <- 10
  for(i in 10:nrow(winmat)) {
    stab <- rbind(stab,data.frame(rep = i,prop = apply(matrix(winmat[1:i,]),2,sum)/sum(apply(matrix(winmat[1:i,]),2,sum)),model = paste0('model(',1:M,')')))
  }

  #set plots to NA if no plots are requested
  p <- p1 <- p2 <- NA

  if(plots) {
    #plot stability
    p1 <- ggplot(stab,aes(y=prop,x=rep,colour=model))
    p1 <- p1 + geom_line() #+ ggtitle(paste0(my.ylab,' cumulative proportion of wins (',K,'-fold, ',repeats,' repeats)'))
    p1 <- p1 + theme(legend.position="none") + xlab(NULL) + ylab('prop. wins')
    if(gray) p1 <- p1 + scale_color_grey()

    # make boxplot
    p <- ggplot(cv.pe, aes(x=Model, y=RMSEP,fill=factor(Model)))
    p <- p + geom_boxplot(aes(group = factor(Model))) + geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3))
    p <- p + scale_x_continuous(breaks = seq(1,M), sec.axis = sec_axis(trans ~ ., name = "Number of Wins", breaks = seq(1,M), labels = wins))
    p <- p + ylab(my.ylab)
    p <- p + theme(legend.position="none")
    if(gray) p <- p + scale_fill_grey(start=.3,end=.7)

    #make density plots
    #detect constants
    cons <- by(cv.pe$RMSEP,as.factor(cv.pe$Model),function(x) !as.logical(sum(diff(x))))

    #plot
    p2 <- ggplot(cv.pe,aes(x=RMSEP,fill=factor(Model),cut=factor(Model))) + theme(legend.position="none",axis.ticks.y=element_blank(),axis.text.y = element_blank())
    for(i in which(cons)) p2 <- p2 + geom_vline(xintercept = cv.pe$RMSEP[cv.pe$Model==i][1])
    p2 <- p2 + geom_density(alpha = 0.5)
    p2 <- p2 + scale_y_continuous(name = ' ',breaks = 0, labels = ' ', sec.axis = sec_axis(~ .,name=' ',breaks = 0, labels = ' '))
    p2 <- p2 + coord_flip() + xlab(NULL) + ylab(NULL)
    if(gray) p2 <- p2 + scale_fill_grey() + scale_color_grey()

    #make titles
    titleplot <- ggplot() + geom_point(aes(1,1), colour="white") +
      theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank()) + annotate("text", x=1, y=1,
        label=paste0(my.ylab,'\n (',K,'-fold, ',repeats,' repeats) \nModel ',which.max(wins),' wins.'))

    #plot all
    grid.arrange(p1, titleplot, p, p2, ncol=2, nrow=2, widths=c(5, 2), heights=c(2, 5))
  }

  #make glm lists
  glmlist <- list()
  for(m in 1:M) {
    glmlist[[m]] <- glm(models[[m]], family = glm.family, data = data)
  }

  #make console output
  linelen <- 60
  prop <- wins/sum(wins)
  l3 <- c(paste0('Results for (',K,'-fold, ',repeats,' repeats)\n'))
  l3 <- c(l3, paste0(' Model:',paste0(rep(' ',linelen),collapse=''),'  |   Wins   |    2.5% |    mean |   97.5% |\n'))

  for(i in 1:length(models)) {
    rmsepci <- quantile(cv.pe$RMSEP[cv.pe$Model==i],c(.025,.975))
    rmsepmean <- mean(cv.pe$RMSEP[cv.pe$Model==i])

    modstring <- deparse(models[[i]],width.cutoff = linelen)[1]
    mx <- nchar(modstring)
    if(mx>linelen) mx <- linelen
    modstring <- substr(modstring,1,mx)
    space <- paste0(rep(' ',(linelen+1) - nchar(modstring)),collapse='')
    l3 <- c(l3, c(paste0(sprintf(' [%2d] ',i),modstring,space,'  |  ', sprintf('%4s',as.character(round(prop[i]*100))),'%   |',sprintf('%8.3f |',round(rmsepci[1],3)),sprintf('%8.3f |',round(rmsepmean,3)),sprintf('%8.3f |',round(rmsepci[2],3)),'\n')))
  }

  #make output list
  output <- list(
    models = models,
    glms = glmlist,
    data = data,
    seed = seed,
    preds = list(preds = preds,y = y),
    stab.plot = p1,
    box.plot = p,
    den.plot = p2,
    win.matrix = winmat,
    wins = wins,
    summary = l3,
    RMSEP = cv.pe)
  attr(output,"class") <- 'xval.glm'

  #output to console
  if(showConsoleOutput) cat(output$summary)

  return(invisible(output))

}

#' cross validation function (only used internally)
cross.validate <- function(M, folds, n, K, glm.family, data, y, models, loss) {

    #set RMSEP
    RMSEP <- numeric(M)

    #resample folds
    folds <- sample(folds, replace = FALSE)

    total_pred <- matrix(NA,n,M)

    #loop of models
    for(m in 1:M) {

      preds <- matrix(NA,n,1)

      #loop of folds
      for(k in 1:K) {
        idx <- which(folds==k,arr.ind=TRUE)
        out <- glm(models[[m]], family = glm.family, data = data[-idx,]) # fit model K-1 folds
        preds[idx,1] <- predict(out, newdata = data[idx,], type = "response") # predict on hold-out
      }

      #calculate loss function
      if(is.null(loss)) {
        RMSEP[m] = sqrt(mean((y - preds)^2))
      } else {
        RMSEP[m] = loss(y,preds)
      }

      total_pred[,m] <- preds
    }

    return(list(loss = RMSEP,pred = total_pred,data = y))
}


#' @export
print.xval.glm <- function(x,...) {
  cat(x$summary)
}


