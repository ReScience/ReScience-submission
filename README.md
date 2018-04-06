
### ReScience submission repository

This is the submission repository for the [Re**Science** journal](https://rescience.github.io).


### Reference implementation of: 

> The cell-type specific cortical microcircuit: relating structure and activity in a full-scale spiking network model. (2014)

The source code is written in Python (version 2.7.5) Brian (version 2.0.1) using Numpy (version
1.11.1) Pandas (version 0.19.2) and Matplotlib (version 2.0.0).

* **article** -> This folder contains the accompanying text.
* **code** -> Here is the source code.
* **data** -> This folder is empty. 
* **notebook** -> There is no notebook for this implementation. 


### How to submit an article ?


1. Create a [github](https://github.com) account

2. [Fork](https://help.github.com/articles/fork-a-repo/) the [ReScience submission](https://github.com/ReScience/ReScience-submission) repository

3. Clone this new repository into your desktop environment

   ```
   $ git clone https://github.com/YOUR-USERNAME/ReScience-submission
   ```

4. Create a branch (the branch name should be author names separated with dashes)

   ```
   $ git checkout -b AUTHOR1-AUTHOR2
   ```


5. Add your code & article (see [author guidelines](https://rescience.github.io/write)) and commit your changes:

   ```
   $ git commit -a -m "Some comment"
   ```


6. [Push](https://help.github.com/articles/pushing-to-a-remote/) to github

   ```
   $ git push origin AUTHOR1-AUTHOR2
   ```

7. Issue a [pull request](https://help.github.com/articles/using-pull-requests/) (PR) to Re**Science** with title "Review Request" and insert the following text in the description:

  ```
  **AUTHOR**

  Dear @ReScience/editors,

  I request a review for the reproduction of the following paper:

  * References of the paper holding results you're replicating

  I believe the original results have been faithfully reproduced as explained in the accompanying article.
  ```

8. Assign the PR to an editor from the [editorial board](https://rescience.github.io/board).

9. Answer questions and requests made in the PR conversation page.

