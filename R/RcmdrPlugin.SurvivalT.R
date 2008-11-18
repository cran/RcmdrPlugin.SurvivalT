.packageName <- "RcmdrPlugin.SurvivalT"
# Rcmdr dialogs, and aditional functions for the RcmdrPlugin.SurvivalT package

# last modified: 19 11 2008 by Daniel Leucuta

.First.lib <- function(libname, pkgname){# copy from TeachingDemos
    if (!interactive()) return()
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        closeCommander(ask=FALSE, ask.save=TRUE)
        Commander()
        }
    }

fncLogRankTest <- function(){ # based on twoSampleWilcoxonTest from Rcmdr 
    initializeDialog(title=gettextRcmdr("Log-rank Test"))
    groupBox <- variableListBox(top, TwoLevelFactors(), title=gettextRcmdr("Groups (pick one)"))
    SurvivalTimeBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Survival Time Variable (pick one)"))
    StatusBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Status Variable (pick one)"))
    onOK <- function(){
        group <- getSelection(groupBox)
        if (length(group) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a groups variable."))
            return()
            }
        Status <- getSelection(StatusBox)
        if (length(Status) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a status variable."))
            return()
            }
        SurvivalTime <- getSelection(SurvivalTimeBox)
        if (length(SurvivalTime) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a survival time variable."))
            return()
            }
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        logger("Censoring status:")
        doItAndPrint(paste("table(", paste(.activeDataSet,"$", Status, sep=""),
            ", ", paste(.activeDataSet,"$", group, sep=""), ")", sep=""))
        logger("Survival time description:")
        doItAndPrint(paste("tapply(", paste(.activeDataSet,"$", SurvivalTime, sep=""),
            ", ", paste(.activeDataSet,"$", group, sep=""), ", summary, na.rm=TRUE)", sep=""))
        logger("Test results:")
        doItAndPrint(paste("survdiff(Surv(", SurvivalTime, " , ", Status, ")~", group, 
				    ', data=', .activeDataSet, ")", sep=""))
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="survdiff")
    tkgrid(getFrame(groupBox), getFrame(SurvivalTimeBox), getFrame(StatusBox), sticky="nw")
    groupsLabel(groupsBox=groupBox, columnspan=2)
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }  

listCoxModels <- function(envir=.GlobalEnv, ...) {  # based on listGeneralizedLinearModels from Rcmdr
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
#        function(.x) "coxph" == (class(eval(parse(text=.x), envir=envir))[1]))]
        function(.x) "coxph" == (class(get(.x, envir=envir))[1]))]
    }

# based on lmP from Rcmdr
CoxP <- function() activeModelP() && class(get(ActiveModel()))[1] == 'coxph'
    
fncCoxModel <- function(){   # based on linearModel from Rcmdr
    # add the class coxph to the modelClasses
    xx <- getRcmdr("modelClasses")
    bolCoxphExists = FALSE
    for(ii in 1:length(xx)){if (xx[ii] == "coxph") bolCoxphExists = TRUE}
    if (bolCoxphExists == FALSE) putRcmdr("modelClasses", c(getRcmdr("modelClasses"), "coxph"))

    #logger("before initialize")
    initializeDialog(title=gettextRcmdr("Cox Model"))
    #logger("before active")
    .activeModel <- ActiveModel()
    #logger("before current if")

    currentModel <- if (!is.null(.activeModel))
        class(get(.activeModel, envir=.GlobalEnv))[1] == "coxph"
#        eval(parse(text=paste("class(", .activeModel, ")[1] == 'lm'", sep="")),
#            envir=.GlobalEnv)
        else FALSE
currentModel <- FALSE
   #logger("before if current")
 #   if (currentModel) {
#        currentFields <- formulaFieldsCox(get(.activeModel, envir=.GlobalEnv))
#        currentFields <- formulaFields(eval(parse(text=.activeModel),
#            envir=.GlobalEnv))
#        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#        }
    #logger("before updateMNr")
    UpdateModelNumber()
    modelName <- tclVar(paste("CoxModel.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()
        if (!is.valid.name(modelValue)){
            errorCondition(recall=fncCoxModel, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
        check.empty <- gsub(" ", "", tclvalue(SurvivalTimeVariable))
        if ("" == check.empty) {
            errorCondition(recall=fncCoxModel, message=gettextRcmdr("Survival time variable of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(StatusVariable))
        if ("" == check.empty) {
            errorCondition(recall=fncCoxModel, message=gettextRcmdr("Status variable of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=fncCoxModel, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
            return()
            }
        if (is.element(modelValue, listCoxModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                fncCoxModel()
                return()
                }
            }
        formula <- paste("Surv(", tclvalue(SurvivalTimeVariable), ", ", tclvalue(StatusVariable), ")~ ", tclvalue(rhsVariable), sep="")
        command <- paste("coxph(", formula,
            ", data=", ActiveDataSet(), subset, ")", sep="")
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        logger("The Cox model: ")
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="coxph", model=TRUE)
    tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    modelFormulaCox()
    subsetBox(model=TRUE)
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }

fncCoxPHSchoenfeld <- function(){ # based on summarizeModel from Rcmdr
    .activeModel <- ActiveModel()
    #if (is.null(.activeModel) || !checkMethod("summary", .activeModel)) return()
    logger("Cox proportional hazard assumption based on scaled Schoenfeld residuals correlated with a Kaplan Meier transformation of time:")
    doItAndPrint(paste("cox.zph(", .activeModel, ")", sep=""))
    }

fncCoxPHSchoenfeldGraph <- function(){ # based on summarizeModel from Rcmdr
    .activeModel <- ActiveModel()
    #if (is.null(.activeModel) || !checkMethod("summary", .activeModel)) return()
    logger("Cox proportional hazard assumption based on scaled Schoenfeld residuals correlated with a Kaplan Meier transformation of time:")
    doItAndPrint(paste("ResidualsDFBeta <- residuals(", .activeModel, ", type='dfbeta')", sep=""))
    doItAndPrint(paste("NumberCoeficients <- length(names(coef(", .activeModel, ")))", sep=""))  
		doItAndPrint(paste("NumberRows <- round(NumberCoeficients+0.1/2) ", sep=""))  
    doItAndPrint(paste("if (NumberRows >= 2) {par(mfrow=c(NumberRows,2))}  
      plot(cox.zph(", .activeModel,"))", sep=""))
    }

fncCoxInfluencePlot <- function(){ # based on summarizeModel from Rcmdr
    #.activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("summary", .activeModel)) return()
    logger("Cox regression index plots of dfbetas to identify influential observations:")
    doItAndPrint(paste("ResidualsDFBeta <- residuals(", .activeModel, ", type='dfbeta')", sep=""))
    doItAndPrint(paste("NumberCoeficients <- length(names(coef(", .activeModel, ")))", sep=""))  
		doItAndPrint(paste("NumberRows <- round(NumberCoeficients+0.1/2) ", sep=""))  
    doItAndPrint(paste("if (NumberRows >= 2) {
		    par(mfrow=c(NumberRows,2))
 				for(i in 1:NumberCoeficients) {
	          plot(ResidualsDFBeta[,i], ylab=names(coef(", .activeModel, "))[i])
	          abline(h=0, lty=3)
        } 
			} else {
				plot(ResidualsDFBeta, ylab=names(coef(", .activeModel, "))[1])
	      abline(h=0, lty=3)
			}	
				", sep=""))
    }

formulaFieldsCox <- function(model){ # based on formulaFields from Rcmdr
    #logger("before formula, after formulaFieldsCox entering")
    formula <- as.character(model$call$formula)
    rhs <- formula[3]
    #logger("before data")
    #SurvivalTimeVariable <- formula[2]
    #StatusVariable <- formula[4] 
    data <- as.character(model$call$data)
    #logger("before decompose")
    SurvivalTimeVariable <- (decomposeSurv(model,data=data)[3]$resp[,2])
    StatusVariable <- (decomposeSurv(model,data=data)[3]$resp[,3])
    #logger("before which.subset")
    which.subset <- which("subset" == names(model$call))
    #logger("before subset")
    subset <- if (0 == length(which.subset)) ""
        else as.character(model$call)[[which.subset]]
    #if (glm) {
    #    fam <- as.character(model$call$family)
    #    family <- fam[1]
    #    link <- fam[2]
    #    }
    #else {
    #    family <- NULL
    #    link <- NULL
    #    }
    doItAndPrint(paste("SurvivalTimeVariable ", SurvivalTimeVariable, "StatusVariable=", StatusVariable, "rhs=", rhs, sep=""))  
    list(SurvivalTimeVariable=SurvivalTimeVariable, StatusVariable=StatusVariable, rhs=rhs, data=data, subset=subset, family=family, link=link)  # !!!!!!!!!!!!!!!
    }
    
modelFormulaCox <- defmacro(frame=top, hasLhs=TRUE, expr={  # based on modelFormula from Rcmdr
    checkAddOperator <- function(rhs){
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1) return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
                rhs.chars[1] else rhs.chars[2]
        !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
        }
    .variables <- Variables()
    word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep="")
    variables <- paste(.variables,
        ifelse(is.element(.variables, Factors()), paste("[", gettextRcmdr("factor"), "]", sep=""), ""))
    xBox <- variableListBox(frame, variables, title=gettextRcmdr("Variables (double-click to formula)"))
    onDoubleClick <- if (!hasLhs){
        function(){
            var <- getSelection(xBox)
            if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
            tkfocus(rhsEntry)
            rhs <- tclvalue(rhsVariable)
            rhs.chars <- rev(strsplit(rhs, "")[[1]])
            check.char <- if (length(rhs.chars) > 0){
                if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
                    rhs.chars[1] else rhs.chars[2]
                }
                else ""
            tclvalue(rhsVariable) <- if (rhs == "" ||
                is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                    paste(rhs, var, sep="")
                else paste(rhs, "+", var)
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
            }
        }
    else{
        function(){
            var <- getSelection(xBox)
            if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
            lhs <- tclvalue(SurvivalTimeVariable)
            lhs2 <- tclvalue(StatusVariable)            
            if (lhs == "") tclvalue(SurvivalTimeVariable) <- var
            else 
            {
            if (lhs2 == "") tclvalue(StatusVariable) <- var
            else {
                tkfocus(rhsEntry)
                rhs <- tclvalue(rhsVariable)
                rhs.chars <- rev(strsplit(rhs, "")[[1]])
                check.char <- if (length(rhs.chars) > 0){
                    if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
                        rhs.chars[1] else rhs.chars[2]
                    }
                    else ""
                tclvalue(rhsVariable) <- if (rhs == "" ||
                    is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                        paste(rhs, var, sep="")
                    else paste(rhs, "+", var)
                }
            }    
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
            }
        }
    tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
    onPlus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "+ ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onTimes <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "*", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onColon <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ":", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onSlash <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onIn <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "%in% ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "- ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onPower <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onLeftParen <- function(){
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onRightParen <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    outerOperatorsFrame <- tkframe(frame)
    operatorsFrame <- tkframe(outerOperatorsFrame)
    plusButton <- buttonRcmdr(operatorsFrame, text="+", width="3", command=onPlus)
    timesButton <- buttonRcmdr(operatorsFrame, text="*", width="3", command=onTimes)
    colonButton <- buttonRcmdr(operatorsFrame, text=":", width="3", command=onColon)
    slashButton <- buttonRcmdr(operatorsFrame, text="/", width="3", command=onSlash)
    inButton <- buttonRcmdr(operatorsFrame, text="%in%", width="5", command=onIn)
    minusButton <- buttonRcmdr(operatorsFrame, text="-", width="3", command=onMinus)
    powerButton <- buttonRcmdr(operatorsFrame, text="^", width="3", command=onPower)
    leftParenButton <- buttonRcmdr(operatorsFrame, text="(", width="3", command=onLeftParen)
    rightParenButton <- buttonRcmdr(operatorsFrame, text=")", width="3", command=onRightParen)

    tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, minusButton,
        powerButton, leftParenButton, rightParenButton, sticky="w")
    formulaFrame <- tkframe(frame)
    if (hasLhs){
        tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Model Formula:     "), fg="blue"), operatorsFrame)
        SurvivalTimeVariable <- if (currentModel) tclVar(currentFields$SurvivalTimeVariable) else tclVar("")
        StatusVariable <- if (currentModel) tclVar(currentFields$StatusVariable) else tclVar("")
        rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
        rhsEntry <- ttkentry(formulaFrame, width="50", textvariable=rhsVariable)
        rhsXscroll <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(rhs, ...))
        tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
        lhsEntry <- ttkentry(formulaFrame, width="10", textvariable=SurvivalTimeVariable)
        lhsScroll <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(lhsEntry, ...))
        tkconfigure(lhsEntry, xscrollcommand=function(...) tkset(lhsScroll, ...))
        lhsEntry2 <- ttkentry(formulaFrame, width="10", textvariable=StatusVariable)
        lhsScroll2 <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(lhsEntry2, ...))
        tkconfigure(lhsEntry2, xscrollcommand=function(...) tkset(lhsScroll2, ...))
        tkgrid(labelRcmdr(formulaFrame, text="Time: "), lhsEntry, labelRcmdr(formulaFrame, text=", Status: "), lhsEntry2, labelRcmdr(formulaFrame, text=" ~    "), rhsEntry, sticky="w")
        tkgrid(labelRcmdr(formulaFrame, text=""), lhsScroll, labelRcmdr(formulaFrame, text=""), lhsScroll2, labelRcmdr(formulaFrame, text=""), rhsXscroll, sticky="w")
        tkgrid.configure(lhsScroll, sticky="ew")
        }
    else{
        rhsVariable <- tclVar("")
        rhsEntry <- ttkentry(formulaFrame, width="50", textvariable=rhsVariable)
        rhsXscroll <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(rhs, ...))
        tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
        tkgrid(labelRcmdr(formulaFrame, text="   ~ "), rhsEntry, sticky="w")
        tkgrid(labelRcmdr(formulaFrame, text=""), rhsXscroll, sticky="w")
        }
    tkgrid.configure(rhsXscroll, sticky="ew")
    })
    
fncKaplanMeier <- function(){ # based on twoSampleWilcoxonTest from Rcmdr
    initializeDialog(title=gettextRcmdr("Kaplan Meier graph"))
    groupBox <- variableListBox(top, TwoLevelFactors(), title=gettextRcmdr("Groups (pick one)"))
    SurvivalTimeBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Survival Time Variable (pick one)"))
    StatusBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Status Variable (pick one)"))
    onOK <- function(){
        group <- getSelection(groupBox)
        if (length(group) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a groups variable."))
            return()
            }
        Status <- getSelection(StatusBox)
        if (length(Status) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a status variable."))
            return()
            }
        SurvivalTime <- getSelection(SurvivalTimeBox)
        if (length(SurvivalTime) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a survival time variable."))
            return()
            }
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        GroupsLevels <- levels(factor(paste(.activeDataSet,"$", group, sep="")))
        #for(i in 1:length(Groups))
        logger("Censoring status:")
        doItAndPrint(paste("table(", paste(.activeDataSet,"$", Status, sep=""),
            ", ", paste(.activeDataSet,"$", group, sep=""), ")", sep=""))
        logger("Survival time description:")
        doItAndPrint(paste("tapply(", paste(.activeDataSet,"$", SurvivalTime, sep=""),
            ", ", paste(.activeDataSet,"$", group, sep=""), ", summary, na.rm=TRUE)", sep=""))
        logger("Kaplan Meier graph:")
        model <- paste("survfit(Surv(", SurvivalTime, " , ", Status, ")~", group, 
				    ', data=', .activeDataSet, ")", sep="")
        doItAndPrint(paste("plot(", model, ", xlab='Time ()', ylab='Survival Probability',
				    lty=c(1,2,3,4,5,6), legend.pos = 1, legend.text=levels(factor(", .activeDataSet,"$", group, ")))", sep = ""))
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="survfit")
    tkgrid(getFrame(groupBox), getFrame(SurvivalTimeBox), getFrame(StatusBox), sticky="nw")
    groupsLabel(groupsBox=groupBox, columnspan=2)
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }
    
fncLogLogSurvival <- function(){ # based on twoSampleWilcoxonTest from Rcmdr
    initializeDialog(title=gettextRcmdr("Log Log Survival graph"))
    groupBox <- variableListBox(top, TwoLevelFactors(), title=gettextRcmdr("Groups (pick one)"))
    SurvivalTimeBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Survival Time Variable (pick one)"))
    StatusBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Status Variable (pick one)"))
    onOK <- function(){
        group <- getSelection(groupBox)
        if (length(group) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a groups variable."))
            return()
            }
        Status <- getSelection(StatusBox)
        if (length(Status) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a status variable."))
            return()
            }
        SurvivalTime <- getSelection(SurvivalTimeBox)
        if (length(SurvivalTime) == 0) {
            errorCondition(recall=fncLogRankTest, message=gettextRcmdr("You must select a survival time variable."))
            return()
            }
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        GroupsLevels <- levels(factor(paste(.activeDataSet,"$", group, sep="")))
        #for(i in 1:length(Groups))
        logger("Censoring status:")
        doItAndPrint(paste("table(", paste(.activeDataSet,"$", Status, sep=""),
            ", ", paste(.activeDataSet,"$", group, sep=""), ")", sep=""))
        logger("Survival time description:")
        doItAndPrint(paste("tapply(", paste(.activeDataSet,"$", SurvivalTime, sep=""),
            ", ", paste(.activeDataSet,"$", group, sep=""), ", summary, na.rm=TRUE)", sep=""))
        logger("log(-log(Survival Probability)) graph:")
        model <- paste("survfit(Surv(", SurvivalTime, " , ", Status, ")~", group, 
				    ', data=', .activeDataSet, ")", sep="")
        doItAndPrint(paste("plot(", model, ", xlab='Time ()', ylab='Log(-Log(Survival Probability))',
				    lty=c(1,2,3,4,5,6), legend.pos = 1, legend.text=levels(factor(", .activeDataSet,"$", group, ")), fun='cloglog')", sep = ""))
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="survfit")
    tkgrid(getFrame(groupBox), getFrame(SurvivalTimeBox), getFrame(StatusBox), sticky="nw")
    groupsLabel(groupsBox=groupBox, columnspan=2)
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }
    