
A Simple BSP-based Model to Predict Execution Time in GPU Applications
======
Authors: Marcos Amarís, Daniel Cordeiro, Alfredo Goldman and Raphael Camargo
======

Models are useful to represent abstractions of software and hardware processes. The Bulk Synchronous Parallel (BSP) is a bridging model for parallel computation that allows algorithmic analysis of programs on parallel computers using performance modeling. The main idea of BSP model is the treatment of communication and computation as abstractions of a parallel system. Meanwhile, the use of GPU devices are becoming more widespread and they are currently capable of performing efficient parallel computation for applications that can be decomposed on thousands of simple threads. However, few models for predicting application execution time on GPUs have been proposed.

In this work we present a simple and intuitive BSP-based model for predicting the CUDA application execution times on GPUs. The model is based on the number of computations and memory accesses of the GPU, with additional information on cache usage obtained from profiling. Scalability, divergence, effect of optimizations and differences of architectures are adjusted by a single parameter.

We evaluated our model using two applications and six different boards. We showed by using profile information for a single board, that the model is general enough to predict the execution time of an application with different input sizes and on different boards with the same architecture. Our model predictions were within $0.8$ to $1.2$ times the measured execution times, which are reasonable for such a simple model. These results indicate that the model is good enough to generalize the predictions for different problem sizes and GPU configurations.
