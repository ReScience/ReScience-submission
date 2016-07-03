
### ReScience submission repository

This is the submission repository for the [Re**Science** journal](https://rescience.github.io).

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

 
  **AUTHOR**

  Dear @ReScience/editors,

  I request a review for the reproduction of the following paper:

  * Benincà, E., Huisman, J., Heerkloss, R., Jöhnk, K.D., Branco, P., Van Nes, E.H., Scheffer, M. & Ellner, S.P. (2008) *Chaos in a long-term experiment with a plankton community.* Nature, 451, 822–825.

  I believed the original results have been adequately reproduced as explained in the accompanying [article](https://github.com/opetchey/ReScience-submission/raw/petchey-plebani-pennekamp-2016/article/article.pdf).
  
  The repository lives at [https://github.com/opetchey/ReScience-submission/tree/petchey-plebani-pennekamp-2016](https://github.com/opetchey/ReScience-submission/tree/petchey-plebani-pennekamp-2016)

8. Assign the PR to an editor from the [editorial board](https://rescience.github.io/board).

9. Answer questions and requests made in the PR conversation page.
