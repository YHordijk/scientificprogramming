General documentation for exacorr module

Remarks:
- To use the standalone version of ExaCorr use the exacorr.x executable in connection with exacorr.inp (only need to specify EXACC section)

To Do:
- [ ] Implement EOM (@stanpa)
  - [ ] Implement Davidson for EOM (@gomes.asp)
  - [ ] Move Davidson overlap to caar-branch
- [ ] Create module for parameters

<table>
  <tr>
    <th>Input parameters</th><th>Comment</th>
  </tr>
  <tr>
    <td>.PRINT</td><td>Print level, helpful for debugging, be careful > 10 tensors are printed </td>
  </tr>       
  <tr>
  	<td>.OCCUPIED</td><td>As list of MOs or energy range  </td>
  </tr>
  <tr>
  	<td>.VIRTUAL</td><td>As list of MOs or energy range  </td>
  </tr>
  <tr>
    <td>.OCC_BETA</td><td> list of beta spin MOs or energy range, alpha spin defined by .OCCUPIED, only use in open shell computations  </td>
  </tr>
  <tr>
    <td>.VIR_BETA</td><td>list of beta spin MOs or energy range, alpha spin defined by .VIRTUAL, only use in open shell computations  </td>
  </tr>
  <tr>
    <td>.LSHIFT</td><td>Level shift </td>
  </tr>
  <tr>
    <td>.TALSH_BUFF</td><td>Maximum memory used in TALSH / for density matrix, should be about 80 % of avialable memory</td>
  </tr>
  <tr>
  	<td>.EXATENSOR</td><td>Calculation with ExaTensor, default TALSH</td>
  </tr>
  <tr>
  	<td>.CCDOUBLES</td><td>CCD calculation</td>
  </tr>
  <tr>
    <td>.CC2</td><td>CC2 calculation</td>
  </tr>
  <tr>
    <td>.NOTRIPLES</td><td>keyword to deactivate triples (recommended for ExaTENSOR) </td>
  </tr>
  <tr>
  	<td>.TCONVERG</td><td>Set convergence criterria (CC iterations, Lambda equations) </td>
  </tr>
  <tr>
  	<td>.LAMBDA</td><td>Solve Λ-equations, needs to be activated for properties</td>
  </tr>
  <tr>
  	<td>.EXA_BLOCKSIZE</td><td>Number to branch the tensors, should be smaller than nocc</td>
  </tr>
  <tr>
  	<td>.MOINT_SCHEME</td><td> different schemes for integral transfomation</td>
  </tr>
  <tr>
  	<td> .NCYCLES</td><td>Set the number of CC cycles</td>
  </tr>
  <tr>
    <td> .CHOLESKY</td><td>threshold for cholesky decomposition (.MOINT_SCHEME 42)</td>
  </tr>
</table>
