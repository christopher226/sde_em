SDEのオイラー丸山近似

使い方

int main(void)
{
  struct ah_params_s ah_params;
  ah_params.alpha=1;
  ah_params.theta=1;
  SDE ah(1,&ah_v0,&ah_v1,&ah_params);
}
 

