  ###############################################################
  #参数设置
  # DiscValid$cpgint <- paste(DiscValid$cpg1,DiscValid$cpg2,sep = "*");
  # my.variable.list <- DiscValid$cpgint;
  # 
  # my.in.variable <- c("age","ds_adj1","ds_adj2","ds_adj3","ds_adj4","smk_adj1","smk_adj2","sex_adj_F",
  #                   "stage_adj_F","strata(ctype)");
  # T1 = NULL;T2 = NULL;
  # Time = "OS";
  # Status = "dead";
  # variable.list = my.variable.list;
  # in.variable = my.in.variable;
  # data = alldatatemp;
  # sle = 0.10;
  # sls = 0.11;
  # vif.threshold = 999;
  #参数设置
  
  
  
  # variable.list = DiscValid$cpgint;
  # in.variable = "NULL";
  # T1 = NULL;
  # T2 = NULL;
  # Time = "OS";
  # Status = "dead";
  # data = alldatatemp;
  # sle = 0.05;
  # sls = 0.06;
  # vif.threshold = 500;
  # 
  
  # variable.list = variableList;
  # in.variable = c("age","ds_adj1","ds_adj2","ds_adj3","smk_adj1","smk_adj2","sex_adj_F",
  #                 "stage_adj_F","strata(ctype)","cg07191189",
  #                 "cg01850110","cg23798714","cg05117823","cg16711866",
  #                 "cg00097800","cg01412762","cg21211213","cg03037684",
  #                 "cg05002580","cg20291674","cg23057936","cg24476103",
  #                 "cg26585560","cg14263093","cg16236779","cg24463290",
  #                 "cg26461637","cg07631435","cg10210919","cg05490591",
  #                 "cg05889881");
  # T1 = NULL;
  # T2 = NULL;
  # Time = "OS";
  # Status = "dead";
  # data = alldatatemp;
  # sle = 0.10;
  # sls = 0.11;
  # vif.threshold = 1999;
  ###############################################################
  

  
#本函数主要用于逐步回归思想,建立多因素的交互作用模型
#LR:模型的似然比检验
#似然比检验两个模型为:旧变量+新变量主效应 VS 旧变量+新变量主效应与交互效应
#VIF:方差膨胀因子,检验模型复共线性
#本函数增加一条限制:多因素模型系数方向与单因素模型系数方向不一致则不纳入模型
#纳入标准略小于剔除标准,否则模型可能为死循环


CoxInteraStepwise <- function(
  variable.list = DiscValid$cpgint,
  in.variable = "NULL",
  T1 = NULL,
  T2 = NULL,
  Time = "OS",
  Status = "dead",
  data = alldatatemp,
  sle = 0.10,
  sls = 0.11,
  vif.threshold = 999
) {
  
  #三种模型：
  #temp.model只含有旧变量对
  #model_noint加入变量的主效应
  #model加入新变量的主效应和交互效应
  #anova似然比检验计算model与model_noint的差别
  library(car)
  library(survival)
  library(stringr)
  univar.pvalue <- NULL
  temp.model <- NULL
  if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time, 
                                            ", ", Status, ") ~ ", paste(in.variable, collapse = "+"), 
                                            sep = "")), data = data, method = "efron")
  }else if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time, 
                                            ", ", Status, ") ~ ", paste(in.variable, collapse = "+"), 
                                            sep = "")), data = data, method = "efron")
  }else if (is.null(Time)) {
    initial.model <- coxph(as.formula(paste("Surv(", T1, 
                                            ", ", T2, ", ", Status, ") ~ ", paste(in.variable, 
                                                                                  collapse = "+"), sep = "")), data = data, method = "efron")
  }else if (is.null(Time)) {
    initial.model <- coxph(as.formula(paste("Surv(", T1, 
                                            ", ", T2, ", ", Status, ") ~ ", paste(in.variable, 
                                                                                  collapse = "+"), sep = "")), data = data, method = "efron")
  }

  if (is.null(initial.model$coefficients)) {
    for (i in 1:length(variable.list)) {
      if (is.null(T2)) {
        uni.model <- coxph(as.formula(paste("Surv(", 
                                            Time, ", ", Status, ") ~ ", variable.list[i], 
                                            sep = "")), data = data, method = "efron")
      }
      if (is.null(Time)) {
        uni.model <- coxph(as.formula(paste("Surv(", 
                                            T1, ", ", T2, ", ", Status, ") ~ ", variable.list[i], 
                                            sep = "")), data = data, method = "efron")
      }
      univar.pvalue[i] <- summary(uni.model)$coefficients[5]
    }
    variable.list1 <- variable.list[univar.pvalue <= 0.9 & 
                                      !is.na(univar.pvalue)]
    univar.pvalue1 <- univar.pvalue[univar.pvalue <= 0.9 & 
                                      !is.na(univar.pvalue)]
    uni.x <- variable.list1[which.min(univar.pvalue1)]
    if (length(uni.x) > 0) {
      if (is.null(T2)) {
        formula <- as.formula(paste("Surv(", Time, ", ", 
                                    Status, ") ~ ", uni.x, sep = ""))
        temp.model <- coxph(formula, data = data, method = "efron")
        if (length(temp.model$coefficients) > 1) {
          print(vif(glm(paste(Status, paste(names(temp.model$coefficients), 
                                            collapse = "+"), sep = "~"), data = data, 
                        family = binomial(link = "logit"))))
        }
      }
      if (is.null(Time)) {
        formula <- as.formula(paste("Surv(", T1, ", ", 
                                    T2, ", ", Status, ") ~ ", uni.x, sep = ""))
        temp.model <- coxph(formula, data = data, method = "efron")
      }
      cat("# --------------------------------------------------------------------------------------------------\n")
      cat("# Initial Model:\n")
      print(summary(temp.model))
    }
  }else if (!is.null(initial.model$coefficients)) {
    temp.model <- initial.model
    cat("# --------------------------------------------------------------------------------------------------\n")
    cat("# Initial Model:\n")
    print(summary(temp.model))
  }#已经运行
  if (length(temp.model$coefficients) > 1) {
    cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
    cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
    print(vif(glm(paste(Status, paste(names(temp.model$coefficients), 
                                      collapse = "+"), sep = "~"), data = data, family = binomial(link = "logit"))))
  }#已经运行，计算方差膨胀因子，探索复共线性

  i <- 0
  break.rule <- TRUE
  while (break.rule) {
    i <- i + 1
    #为防止交互作用项左右顺序调换相反
    if (i == 1) {
      modelvar <- rownames(summary(temp.model)$coefficients)
      variable.list2 <- setdiff(variable.list, gsub(":","*",modelvar,fixed = TRUE))
      variable.list2_int <- c(gsub("*",":",variable.list2,fixed = TRUE))
    }else {
      modelvar <- rownames(summary(temp.model)$coefficients)
      modelvar2 <- stringr::str_split(str_subset(modelvar,"[:]"),":")
      for (ii in 1:length(modelvar2)){
        modelvar2[[ii]] <- paste(modelvar2[[ii]][2],modelvar2[[ii]][1],sep = ":")
      }
      modelvar <- c(modelvar,unlist(modelvar2))
      variable.list2 <- setdiff(variable.list, c(gsub(":","*",modelvar,fixed = TRUE)))
      variable.list2_int <- c(gsub("*",":",variable.list2,fixed = TRUE))
      out.x <- NULL
      out.x.1 <- NULL
    }#注意out.x,已经运行,variable.list2在表达式中不包含的变量
    # variable.list2为*形式
    # variable.list2_int为：形式
    
    #将variable.list2_int交互项左右交换防止因为左右项相反而识别不到
    variable.list2_int_inverse <- stringr::str_split(variable.list2_int,":")
    for (ii in 1:length(variable.list2_int_inverse)){
      variable.list2_int_inverse[[ii]] <- paste(variable.list2_int_inverse[[ii]][2],variable.list2_int_inverse[[ii]][1],sep = ":")
    }
    variable.list2_int_inverse <- unlist(variable.list2_int_inverse)
    

    
    if (length(variable.list2) != 0) {
      anova.pvalue <- NULL
      mv.pvalue <- NULL
      vif.value <- NULL
      wald.betaint <- NULL
      wald.betaintsingle <- NULL
      wald.betamain1 <- NULL
      wald.betamain2 <- NULL
      wald.betamain1single <- NULL
      wald.betamain2single <- NULL
      
      # 修改程序，做只有主效应与同时有主效应和交互的似然比检验
      # 结果衡量交互效应大小对于模型的贡献
      for (k in 1:length(variable.list2)) {
        model <- update(temp.model, as.formula(paste(". ~ . + ", 
                                                     variable.list2[k], sep = "")))
        model_noint <- update(model, as.formula(paste(". ~ . - ", 
                                                      variable.list2_int[k], sep = "")))
        if (length(model$coefficients) > 1) {
          if (sum(is.na(model$coefficients)) != 0) {
            anova.pvalue[k] <- 1
            mv.pvalue[k] <- 1
            vif.value[k] <- 999
            wald.betaint[k] <- NULL
            wald.betaintsingle[k] <- NULL
            wald.betamain1[k] <- NULL
            wald.betamain2[k] <- NULL
            wald.betamain1single[k] <- NULL
            wald.betamain2single[k] <- NULL
          }else {
            anova.pvalue[k] <- anova(model_noint, model)[2, 
                                                         "P(>|Chi|)"]
            mv.pvalue[k] <- summary(model)$coefficients[nrow(summary(model)$coefficients), 
                                                        "Pr(>|z|)"]
            model.vif <- vif(glm(as.formula(paste(Status, 
                                                  paste(names(model$coefficients), collapse = "+"), 
                                                  sep = "~")), data = data, family = binomial(link = "logit")))
            vif.value[k] <- model.vif[length(model.vif)]
            wald.betamain1[k] <- summary(model)$coefficients[,"coef"][rownames(summary(model)$coefficients) %in%  
                                                                        unlist(stringr::str_split(variable.list2_int[k],":"))[1]]
            
            wald.betamain2[k] <- summary(model)$coefficients[,"coef"][rownames(summary(model)$coefficients) %in%  
                                                                        unlist(stringr::str_split(variable.list2_int[k],":"))[2]]
            
            wald.betaint[k] <- summary(model)$coefficients[,"coef"][rownames(summary(model)$coefficients) %in%  
                                                                     c(variable.list2_int[k],variable.list2_int_inverse[k])]
            m <- coxph(as.formula(paste("Surv(", Time, 
                                        ", ", Status, ") ~ ", paste(c(in.variable,variable.list2[k]), collapse = "+"), 
                                        sep = "")), data = data, method = "efron")
            
            wald.betamain1single[k] <- summary(m)$coefficients[,"coef"][rownames(summary(m)$coefficients) %in%  
                                                                        unlist(stringr::str_split(variable.list2_int[k],":"))[1]]
            
            wald.betamain2single[k] <- summary(m)$coefficients[,"coef"][rownames(summary(m)$coefficients) %in%  
                                                                        unlist(stringr::str_split(variable.list2_int[k],":"))[1]]
            wald.betaintsingle[k] <- summary(m)$coefficients[, 
                                                             "coef"][rownames(summary(m)$coefficients) %in% 
                                                                       c(variable.list2_int[k],variable.list2_int_inverse[k])]
          }
        }
      }#for结束
      

      #for 循环结束，一个一个变量放入模型检验,评价模型外变量重要性
      #计算anova.pvalue（似然比模型P值），mv.pvalue（交互项p值），model.vif（模型所有VIF），vif.value（交互项VIF）
      
      #去掉交互项P值大于0.9的交互变量名，anova.pvalue（模型P值），mv.pvalue（交互项p值）
      #选择标准为anova值小于sle
      #如果多因素模型与单因素模型系数不统一则不考虑纳入
      variable.list2.1 <- variable.list2[mv.pvalue <= 0.9 & 
                                           !is.na(mv.pvalue) & vif.value <= vif.threshold &  
                                           wald.betaint*wald.betaintsingle > 0 &
                                           wald.betamain1*wald.betamain1single > 0 &
                                           wald.betamain2*wald.betamain2single > 0]
      anova.pvalue2 <- anova.pvalue[mv.pvalue <= 0.9 & 
                                      !is.na(mv.pvalue) & vif.value <= vif.threshold &  
                                      wald.betaint*wald.betaintsingle > 0 &
                                      wald.betamain1*wald.betamain1single > 0 &
                                      wald.betamain2*wald.betamain2single > 0]
      mv.pvalue2 <- mv.pvalue[mv.pvalue <= 0.9 & !is.na(mv.pvalue) & 
                                vif.value <= vif.threshold &  
                              wald.betaint*wald.betaintsingle > 0 &
                                wald.betamain1*wald.betamain1single > 0 &
                                wald.betamain2*wald.betamain2single > 0]
      if (length(variable.list2.1) == 0){
        enter.x <- NULL
      }else{
        #进入的变量取anova.pvalue（模型P值）最小的变量对(名+wald.p)，且小于sle
        enter.x <- variable.list2.1[anova.pvalue2 == min(anova.pvalue2, 
                                                         na.rm = TRUE) & anova.pvalue2 <= sle]
        wald.p <- mv.pvalue2[anova.pvalue2 == min(anova.pvalue2, 
                                                  na.rm = TRUE) & anova.pvalue2 <= sle]
      }

      
      if (length(setdiff(enter.x, NA)) != 0) {
        if (length(enter.x) > 1) {
          enter.x <- enter.x[which.min(wald.p)]
        }
        cat("# --------------------------------------------------------------------------------------------------", 
            "\n")
        cat(paste("### iter num = ", i, ", Forward Selection by LR Test: ", 
                  "+ ", enter.x, sep = ""), "\n")
        temp.model <- update(temp.model, as.formula(paste(". ~ . + ", 
                                                          enter.x, sep = "")))
        print(summary(temp.model))
        cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
        cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
        if (length(temp.model$coefficients) > 1) {
          print(vif(glm(paste(Status, paste(names(temp.model$coefficients), 
                                            collapse = "+"), sep = "~"), data = data, 
                        family = binomial(link = "logit"))))
        }
      }
      
    }else {
      enter.x <- NULL
    }
    #if (length(variable.list2) != 0)条件语句结束
    #加入一个变量对后重新计算VIF,选取anova小的，anova低于sle的;变量筛选过程,enter.x小于入选标准
    
    # enter.x为*号形式
    if (length(enter.x) != 0 ){
      enter.x_main_int <- unlist(c(strsplit(enter.x, "*", fixed=TRUE),gsub("*",":",enter.x, fixed=TRUE)))
    }else{
      enter.x_main_int <- NULL
    }
    

    # 剔除变量开始
    # 选择模型中的交互变量对（除刚纳入的外），一对一对进行检验，剔除
    # variable.list3变换为*号形式，variable.list3为待检验是否剔除变量列表（除新纳入的外）
    if (i == 1 & length(enter.x) == 0) {
      cat("# ==================================================================================================", 
          "\n")
      cat(paste("*** Stepwise Final Model (in.lr.test: sle = ", 
                sle, "; variable selection restrict in vif = ", 
                vif.threshold, "):", sep = ""), "\n")
      print(summary(temp.model))
      break
    }else {
      #modelvar2 <- rownames(summary(temp.model)$coefficients)
      variable.list3 <- setdiff(rownames(summary(temp.model)$coefficients), 
                                c(enter.x_main_int, in.variable))
      
      variable.list3 <- c(gsub(":","*",variable.list3,fixed = TRUE))
      variable.list3 <- stringr::str_subset(variable.list3,"[*]")
      variable.list3_int <- c(gsub("*",":",variable.list3,fixed = TRUE))
      
      #剔除变量后计算anova,即似然比检验
      if (length(variable.list3) != 0) {
        anova.pvalue <- NULL
        vif.value <- NULL
        wald.betaint <- NULL
        wald.betaintsingle <- NULL
        for (k in 1:length(variable.list3)) {
          # model <- update(temp.model, as.formula(paste(". ~ . - ", 
          #                                              variable.list3[k], sep = "")))
          model_noint <- update(temp.model, as.formula(paste(". ~ . - ", 
                                                             variable.list3_int[k], sep = "")))
          anova.pvalue[k] <- anova(model_noint, temp.model)[2, 
                                                            "P(>|Chi|)"]

        }
        
        # out.x为*号形式
        # out.x_main_int为主效应+交互“：”形式
        # out.x_int为交互“：”形式
        
        out.x <- variable.list3[(anova.pvalue == max(anova.pvalue, 
                                                     na.rm = TRUE) & anova.pvalue > sls)]
        
        # 如果有多个out.x，取最大的交互作用P值
        out.x <- setdiff(out.x, NA)
        if (length(out.x) != 0) {
          if (length(out.x) > 1) {
            out.x.1 <- out.x
            for (j in 1:length(out.x)) {
              out.x[j] <- out.x.1[(length(out.x) - j + 
                                     1)]
            }
            wald.p <- rep(NA, length(out.x))
            for (j in 1:length(out.x)) {
              wald.p[j] <- summary(temp.model)$coefficients[, 
                                                            "Pr(>|z|)"][rownames(summary(temp.model)$coefficients) == 
                                                                          out.x[j]]
            }
            out.x <- out.x[which.max(wald.p)]
          }
          
          temp.model <- update(temp.model, as.formula(paste(". ~ . - ", 
                                                            out.x, sep = "")))
          summary(temp.model)
          if (length(temp.model$coefficients) > 1) {
            vif(glm(paste(Status, paste(names(temp.model$coefficients), 
                                        collapse = "+"), sep = "~"), data = data, 
                    family = binomial(link = "logit")))
          }
        }
      } else {
        out.x <- NULL
      }#out.x=NULL
    }
    
    #无变量可纳入无变量可剔除，终止程序
    if ((length(enter.x) + length(out.x)) == 0) {
      final.model <- temp.model
      cat("# ==================================================================================================", 
          "\n")
      cat(paste("*** Stepwise Final Model (in.lr.test: sle = ", 
                sle, "; out.lr.test: sls = ", sls, "; variable selection restrict in vif = ", 
                vif.threshold, "):", sep = ""), "\n")
      print(summary(final.model))
      cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
      cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
      if (length(final.model$coefficients) > 1) {
        print(vif(glm(paste(Status, paste(names(final.model$coefficients), 
                                          collapse = "+"), sep = "~"), data = data, family = binomial(link = "logit"))))
      }
      break.rule <- FALSE
    }
    enter.x <- NULL
  } #while循环结束
}  




# CoxInteraStepwise(
#   variable.list = DiscValid$cpgint,
#   in.variable = c("age","ds_adj1","ds_adj2","ds_adj3","ds_adj4","smk_adj1","smk_adj2","sex_adj_F",
#                      "stage_adj_F","strata(ctype)"),
#   T1 = NULL,
#   T2 = NULL,
#   Time = "OS",
#   Status = "dead",
#   data = alldatatemp,
#   sle = 0.10,
#   sls = 0.11,
#   vif.threshold = 999
# )

# rm(variable.list,in.variable,T1,T2,Time,Status,data,sle,sls,vif.threshold,vif.value,i,break.rule,
#    univar.pvalue,temp.model,initial.model,variable.list,uni.model,univar.pvalue,variable.list1,uni.x,
#    modelvar,modelvar2,variable.list2,variable.list2_int,out.x,out.x.1,variable.list2_int_inverse,
#    anova.pvalue,mv.pvalue,vif.value,wald.betaint,wald.betaintsingle,model,model_noint,anova.pvalue,
#    mv.pvalue,vif.value,anova.pvalue,mv.pvalue,model.vif,m,wald.betaintsingle,variable.list2.1,
#    anova.pvalue2,mv.pvalue2,enter.x,wald.p,enter.x_main_int,variable.list3,variable.list3_int,
#    anova.pvalue,vif.value,wald.betaint,wald.betaintsingle,anova.pvalue,out.x
# )
