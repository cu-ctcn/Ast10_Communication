# Calculating VIP score using sklearn PLSR model output
import numpy as np
def vip(model):
  print("from: https://github.com/scikit-learn/scikit-learn/issues/7050")
  t = model.x_scores_
  w = model.x_weights_
  q = model.y_loadings_
  p, h = w.shape
  vips = np.zeros((p,))
  s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
  total_s = np.sum(s)
  print(q)
  print("ss: ", total_s)
  print(s)
  print("h:",h)
  for i in range(p):
      weight = np.array([ (w[i,j] / np.linalg.norm(w[:,j]))**2 for j in range(h) ])
      vips[i] = np.sqrt(p*(s.T @ weight)/total_s)
      print("weight: ",weight)
      print("w: ", w[i,:])
      print((s.T @ weight)/total_s)
  return vips
