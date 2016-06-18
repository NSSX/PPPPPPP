/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_lstadd_back.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/20 11:50:55 by igomez            #+#    #+#             */
/*   Updated: 2014/11/20 12:20:06 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

void	ft_lstadd_back(t_list **lst, t_list *n)
{
	t_list	*temp;

	if (!n || !lst)
		return ;
	if (!(*lst))
		*lst = n;
	else
	{
		temp = *lst;
		while (temp->next)
			temp = temp->next;
		temp->next = n;
	}
}
