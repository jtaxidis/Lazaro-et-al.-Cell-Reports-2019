function save_eps(name)

set(gcf,'renderer','painters');
print(gcf,'-depsc2',[name,'.eps'])
