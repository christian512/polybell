import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
times = [195, 14.5]
x = [0,1]
xlabels = ["without symmetries", "with symmetries"]
mybars = plt.bar(x, times, 0.6)
plt.xticks(x, xlabels)
plt.yticks([],[])
# get rid of the frame
for spine in plt.gca().spines.values():
    spine.set_visible(False)
#plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')

# direct label each bar with Y axis values
for bari in mybars:
    height = bari.get_height()
    plt.gca().text(bari.get_x() + bari.get_width()/2, bari.get_height()-12, str(int(height)),
                 ha='center', color='white', fontsize=15)


plt.savefig('calc_times.png', dpi=300)