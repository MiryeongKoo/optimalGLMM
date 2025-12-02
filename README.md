# optimalLMM

## This package supports the automatic exploration of optimal linear mixed models (LMMs) for hierarchical data. 

This algorithm employs a forward selection approach, which starts with a model containing no predictors and then adds predictors one at a time until a complex model with the newly added predictor is not significantly different from a simpler model. Unlike traditional forward selection, which adds the predictor that yields the greatest additional improvement to the model fit, the algorithm selects the predictor based on its statistical importance in prediction, as determined by machine learning methods (here, mixed-effects random forest; Hajjem et al., 2014). The algorithm searches for the best-fitting model through three primary steps: preparation, predictor selection, and model selection. Since exploring the best-fitting LMM is complicated, two variations of the algorithm are developed to reflect the complexity of the optimal LMM selection. 

## Key functions
With a commonly applied preparation step, **optimalLMM** sequentially performs predictor and model selection steps, whereas **optimalLMM2** alternates between these steps. More detailed explanations are provided in the next sections. 

**filter_small_clusters** fits a null mixed-effects model (a random intercept only model without any predictors) using lmer and calculates the intraclass correlation coefficient (ICC) which is the ratio of between-cluster variance to total variance.

**calculate_ICC** removes all observations belonging to clusters that contain fewer than a specified minimum number of subjects.

## Reference
A brief description and example analysis can be found in Koo & Zhang (2025)'s paper, titled [Automatic search algorithm for optimal generalized linear mixed models (GLMMs)](https://aclanthology.org/2025.aimecon-main.38/) in Proceedings of the Artificial Intelligence in Measurement and Education Conference (AIME-Con).
